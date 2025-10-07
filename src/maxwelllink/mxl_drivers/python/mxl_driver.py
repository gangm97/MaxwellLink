#!/usr/bin/env python3

from __future__ import annotations
import argparse
import subprocess
import shlex
import time
import shutil
import os
import socket, json
import numpy as np
import struct

try:
    from .models import __drivers__
    from .models.dummy_model import DummyModel
except ImportError:
    from models import __drivers__
    from models.dummy_model import DummyModel

description = """
A Python driver connecting to MaxwellLink, receiving E-field data and returning
the source amplitude vector for a quantum dynamics model.
"""


# helper function to determine whether this processor is the MPI master using mpi4py
def am_master():
    try:
        from mpi4py import MPI as _MPI

        _COMM = _MPI.COMM_WORLD
        _RANK = _COMM.Get_rank()
    except Exception:
        _COMM = None
        _RANK = 0
    return _RANK == 0


_INT32 = struct.Struct("<i")
_FLOAT64 = struct.Struct("<d")
# numpy dtypes on the wire
DT_FLOAT = np.float64
DT_INT = np.int32

# header width (ASCII, space-padded)
HEADER_LEN = 12

# Message codes similar to i-PI's socket interface
STATUS = b"STATUS"
READY = b"READY"
HAVEDATA = b"HAVEDATA"
NEEDINIT = b"NEEDINIT"
INIT = b"INIT"
POSDATA = b"POSDATA"
GETFORCE = b"GETFORCE"
FORCEREADY = b"FORCEREADY"
STOP = b"STOP"
BYE = b"BYE"

# EM aliases for readability (same wire format)
FIELDDATA = POSDATA
GETSOURCE = GETFORCE
SOURCEREADY = FORCEREADY


class SocketClosed(OSError):
    pass


def _pad12(msg: bytes) -> bytes:
    if len(msg) > HEADER_LEN:
        raise ValueError("Header too long")
    return msg.ljust(HEADER_LEN, b" ")


def send_msg(sock: socket.socket, msg: bytes) -> None:
    """Send a 12-byte ASCII header (space-padded)."""
    sock.sendall(_pad12(msg))


def _recvall(sock: socket.socket, n: int) -> bytes:
    """Read exactly n bytes or raise SocketClosed."""
    buf = bytearray()
    while len(buf) < n:
        chunk = sock.recv(n - len(buf))
        if not chunk:
            raise SocketClosed("Peer closed")
        buf.extend(chunk)
    return bytes(buf)


def recv_msg(sock: socket.socket) -> bytes:
    """Receive 12-byte ASCII header."""
    hdr = _recvall(sock, HEADER_LEN)
    return hdr.rstrip()


def send_array(sock: socket.socket, arr, dtype) -> None:
    a = np.asarray(arr, dtype=dtype, order="C")
    sock.sendall(memoryview(a).cast("B"))


def recv_array(sock: socket.socket, shape, dtype):
    out = np.empty(shape, dtype=dtype, order="C")
    mv = memoryview(out).cast("B")
    need = mv.nbytes
    got = 0
    while got < need:
        r = sock.recv_into(mv[got:], need - got)
        if r == 0:
            raise SocketClosed("Peer closed")
        got += r
    return out


def send_int(sock: socket.socket, x: int) -> None:
    sock.sendall(_INT32.pack(int(x)))


def recv_int(sock: socket.socket) -> int:
    buf = bytearray(_INT32.size)
    mv = memoryview(buf)
    got = 0
    while got < _INT32.size:
        r = sock.recv_into(mv[got:], _INT32.size - got)
        if r == 0:
            raise SocketClosed("Peer closed")
        got += r
    return _INT32.unpack(buf)[0]


def send_bytes(sock: socket.socket, b: bytes) -> None:
    send_int(sock, len(b))
    if len(b):
        sock.sendall(b)


def recv_bytes(sock: socket.socket) -> bytes:
    n = recv_int(sock)
    return _recvall(sock, n) if n else b""


def recv_posdata(sock: socket.socket):
    """Read POSDATA/FIELDDATA block."""
    cell = recv_array(sock, (3, 3), DT_FLOAT).T.copy()
    icell = recv_array(sock, (3, 3), DT_FLOAT).T.copy()
    nat = recv_int(sock)
    xyz = recv_array(sock, (nat, 3), DT_FLOAT)
    return cell, icell, xyz


def send_force_ready(
    sock: socket.socket,
    energy_ha: float,
    forces_Nx3_ha_per_bohr,
    virial_3x3_ha,
    more: bytes = b"",
):
    """FORCEREADY/SOURCEREADY: energy (1), natoms (1), forces (Nx3), virial (3x3), extras (len + bytes)."""
    send_msg(sock, FORCEREADY)
    send_array(sock, np.array([energy_ha], dtype=DT_FLOAT), DT_FLOAT)
    forces = np.asarray(forces_Nx3_ha_per_bohr, dtype=DT_FLOAT)
    assert forces.ndim == 2 and forces.shape[1] == 3
    send_int(sock, forces.shape[0])
    send_array(sock, forces, DT_FLOAT)
    send_array(sock, np.asarray(virial_3x3_ha, dtype=DT_FLOAT).T, DT_FLOAT)
    send_bytes(sock, more)


# the above functions can be also obtained from maxwelllink.sockets,
# but we copy them here to avoid circular imports.


def read_value(s):
    """Attempt to parse a string to int or float; fallback to string."""
    s = s.strip()
    for cast in (int, float):
        try:
            return cast(s)
        except ValueError:
            continue
    if s.lower() == "false":
        return False
    if s.lower() == "true":
        return True
    return s


def read_args_kwargs(input_str):
    """
    Parses a string into positional arguments and keyword arguments.

    Args:
        input_str (str): The input string containing comma-separated values and key-value pairs.

    Returns:
        tuple: A tuple containing a list of positional arguments and a dictionary of keyword arguments.
    """
    args = []
    kwargs = {}
    tokens = input_str.split(",")
    for token in tokens:
        token = token.strip()
        if "=" in token:
            key, value = token.split("=", 1)
            kwargs[key.strip()] = read_value(value)
        elif len(token) > 0:
            args.append(read_value(token))
    return args, kwargs


def run_driver(
    unix=False,
    address="localhost",
    port: int = 31415,
    timeout: float = 600.0,
    driver=DummyModel(),
    sockets_prefix="/tmp/socketmxl_",
):
    """
    Run the driver to communicate with the MaxwellLink simulation via sockets.
    """
    if unix:
        sock = socket.socket(socket.AF_UNIX)
        sock.connect(sockets_prefix + address)
    else:
        sock = socket.socket(socket.AF_INET)
        # NEW: keepalive + nodelay for low-latency, long-lived tiny messages
        try:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
            sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        except (OSError, AttributeError):
            pass

        sock.connect((address, port))
        sock.settimeout(timeout)

    initialized = False
    have_result = False
    pending_amp = None
    additional_data = {}
    dt_au = 0.0
    molid = None

    while True:
        try:
            msg = recv_msg(sock)
        except Exception:
            # Treat EOF/timeouts during normal shutdown as clean exit
            break

        if msg == STATUS:
            # Server is polling; we must reply with our state.
            if not initialized:
                send_msg(sock, NEEDINIT)
            elif have_result:
                send_msg(sock, HAVEDATA)
            else:
                send_msg(sock, READY)

        elif msg == INIT:
            # Server sends INIT after we replied NEEDINIT
            molid = recv_int(sock)
            init_json = json.loads(recv_bytes(sock).decode("utf-8") or "{}")
            dt_au = float(init_json.get("dt_au", 0.0))
            print("[initialization] Time step in atomic units:", dt_au)
            print("[initialization] Assigned a molecular ID:", molid)

            driver.initialize(dt_au, molid)

            initialized = True
            print("[initialization] Finished initialization for molecular ID:", molid)

        elif msg == FIELDDATA or msg == b"POSDATA":
            # One step of data from server: treat "positions" as the E-field vector in a.u.
            # This is to mirror i-pi's existing socket interface.
            cell, icell, xyz = recv_posdata(sock)
            # effective [Ex, Ey, Ez] (a.u.) for this molecule
            E = xyz[0]
            # Stage the step (no commit)
            driver.stage_step(E)

            have_result = True

        elif msg == GETSOURCE or msg == b"GETFORCE":
            # Server asks us to return the result for this step
            if not driver.have_result():
                # it means the driver code was terminated during driver.propagate() and driver.calc_amp_vector()
                # one way is to be defensive: return zero if we somehow got here without a computed result
                pending_amp = np.zeros(3, float)
            else:
                pending_amp = driver.commit_step()
                additional_data = driver.append_additional_data()
            send_force_ready(
                sock,
                energy_ha=0.0,
                forces_Nx3_ha_per_bohr=pending_amp.reshape(1, 3),
                virial_3x3_ha=np.zeros((3, 3)),
                more=json.dumps(
                    additional_data,
                    ensure_ascii=False,
                    separators=(",", ":"),
                    sort_keys=True,
                ).encode("utf-8"),
            )
            have_result = False
            pending_amp = None

        elif msg == STOP:
            # Acknowledge and leave gracefully
            try:
                send_msg(sock, BYE)
            finally:
                print("Received STOP, exiting")
                break
        else:
            raise RuntimeError(f"Unexpected header: {msg!r}")


def mxl_driver_main():
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-u",
        "--unix",
        action="store_true",
        default=False,
        help="Use a UNIX domain socket.",
    )
    parser.add_argument(
        "-a",
        "--address",
        type=str,
        default="localhost",
        help="Host name (for INET sockets) or name of the UNIX domain socket to connect to.",
    )
    parser.add_argument(
        "-S",
        "--sockets_prefix",
        type=str,
        default="/tmp/socketmxl_",
        help="Prefix used for the unix domain sockets. Ignored when using TCP/IP sockets.",
    )
    parser.add_argument(
        "-p",
        "--port",
        type=int,
        default=31415,
        help="TCP/IP port number. Ignored when using UNIX domain sockets.",
    )
    parser.add_argument(
        "-m",
        "--model",
        type=str,
        default="dummy",
        choices=list(__drivers__.keys()),
        help="""Type of molecular / material model for computing dipole moments under EM field.
        """,
    )
    parser.add_argument(
        "-o",
        "--param",
        type=str,
        default="",
        help="""Parameters required to run the driver. Comma-separated list of values
        """,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose output.",
    )

    args = parser.parse_args()

    driver_args, driver_kwargs = read_args_kwargs(args.param)

    if args.model in __drivers__:
        try:
            d_f = __drivers__[args.model](
                *driver_args, verbose=args.verbose, **driver_kwargs
            )
        except ImportError:
            # specific errors have already been triggered
            raise
        except Exception as err:
            print(f"Error setting up molecular dynamics model {args.model}")
            print(__drivers__[args.model].__doc__)
            print("Error trace: ")
            raise err
    elif args.model == "dummy":
        d_f = DummyModel(verbose=args.verbose)
    else:
        raise ValueError("Unsupported driver model ", args.model)

    run_driver(
        unix=args.unix,
        address=args.address,
        port=args.port,
        driver=d_f,
        sockets_prefix=args.sockets_prefix,
    )


def _clean_env_for_subprocess():
    env = os.environ.copy()
    # Nuke anything that makes a child think it's an MPI rank
    prefixes = (
        "PMI_",
        "PMIX_",
        "OMPI_",
        "MPI_",
        "MPICH_",
        "I_MPI_",
        "HYDRA_",
        "SLURM_",
        "FI_",
        "UCX_",
        "PSM2_",
        "PMI",
    )
    for k in list(env.keys()):
        for p in prefixes:
            if k.startswith(p):
                env.pop(k, None)
                break
    # Some MPIs set these exact names without a prefix
    for k in ("PMI_FD", "PMI_PORT", "PMI_ID", "PMI_RANK", "PMI_SIZE"):
        env.pop(k, None)
    return env


def launch_driver(
    command='--model tls --port 31415 --param "omega=0.242, mu12=187, orientation=2, pe_initial=1e-4" --verbose',
    sleep_time=0.5,
):
    """
    Helper function to launch the driver in a local background process.
    Returns the subprocess.Popen object for the launched process.
    """
    if am_master():
        print(f"Launching driver with command: mxl_driver.py {command}")
        # launch the external driver (client)
        driver_argv = shlex.split(shutil.which("mxl_driver.py") + " " + command)
        # Use a fresh, non-blocking subprocess; inherit env/stdio for easy debugging
        proc = subprocess.Popen(driver_argv, env=_clean_env_for_subprocess())
        time.sleep(sleep_time)
        return proc
    else:
        return None


def terminate_driver(proc, timeout=2.0):
    """
    Helper function to terminate the driver process launched by launch_driver.
    """
    if proc is not None and am_master():
        # Give it a moment to shut down naturally after the sim closes the socket
        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.terminate()
            print("Driver did not exit cleanly, sent terminate signal")
            try:
                proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                proc.kill()
                print("Driver did not terminate, sent kill signal")


if __name__ == "__main__":
    mxl_driver_main()
