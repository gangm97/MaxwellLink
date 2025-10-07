"""
Socket layer for MaxwellLink drivers and servers in a single file.
The socket protocol of MaxwellLink is based on i-PI's socket protocol (https://ipi-code.org/).

Public API:
- SocketHub: multi-client server/poller for handling many driver connections with the FDTD engine
- Protocol constants: STATUS, READY, HAVEDATA, NEEDINIT, INIT, ...
- EM aliases: FIELDDATA, GETSOURCE, SOURCEREADY  (map 1:1 to POSDATA/GETFORCE/FORCEREADY)
- Low-level helpers: send_msg, recv_msg, send_array/recv_array, etc.
- Exceptions: SocketClosed
"""

from __future__ import annotations
import socket, struct, json, time, threading, os
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple
import numpy as np

_INT32 = struct.Struct("<i")
_FLOAT64 = struct.Struct("<d")

# Fixed header width (ASCII, space-padded)
HEADER_LEN = 12
# Canonical i-PI message codes
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

# numpy dtypes on the wire (i-PI/ASE use float64 for reals, int32 for counts)
DT_FLOAT = np.float64
DT_INT = np.int32


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


# -------- compound payloads (i-PI compatible) --------


def send_posdata(
    sock: socket.socket, cell_3x3_bohr, invcell_3x3_per_bohr, positions_Nx3_bohr
):
    """POSDATA/FIELDDATA: cell, invcell, natoms, positions (all float64)."""
    assert np.asarray(cell_3x3_bohr).shape == (3, 3)
    assert np.asarray(invcell_3x3_per_bohr).shape == (3, 3)
    pos = np.asarray(positions_Nx3_bohr, dtype=DT_FLOAT)
    assert pos.ndim == 2 and pos.shape[1] == 3
    send_msg(sock, POSDATA)
    send_array(sock, np.asarray(cell_3x3_bohr, dtype=DT_FLOAT).T, DT_FLOAT)
    send_array(sock, np.asarray(invcell_3x3_per_bohr, dtype=DT_FLOAT).T, DT_FLOAT)
    send_int(sock, pos.shape[0])
    send_array(sock, pos, DT_FLOAT)


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


def recv_getforce(sock: socket.socket):
    """On server side after sending GETFORCE/GETSOURCE and receiving FORCEREADY/SOURCEREADY."""
    e = float(recv_array(sock, (1,), DT_FLOAT)[0])
    nat = recv_int(sock)
    frcs = recv_array(sock, (nat, 3), DT_FLOAT)
    vir = recv_array(sock, (3, 3), DT_FLOAT).T.copy()
    extra = recv_bytes(sock)
    return e, frcs, vir, extra


# -------- convenience wrappers for EM (i-PI compatible) --------


def pack_em_fieldata(
    sock: socket.socket, t_au: float, dt_au: float, efield_au_vec3, meta: dict
):
    """Send 'FIELDDATA' as POSDATA with natoms=1 and positions=[Ex,Ey,Ez]."""
    I = np.eye(3, dtype=DT_FLOAT)
    exyz = np.asarray(efield_au_vec3, dtype=DT_FLOAT).reshape(1, 3)
    send_posdata(sock, I, I, exyz)
    # meta/time tags can be sent back in SOURCEREADY's extra blob if needed.


def pack_init(sock: socket.socket, init_dict: dict):
    """Send INIT with JSON payload (molecule_id:int, JSON bytes)."""
    send_msg(sock, INIT)
    molid = int(init_dict.get("molecule_id", 0))
    send_int(sock, molid)
    init_bytes = json.dumps(init_dict).encode("utf-8")
    send_bytes(sock, init_bytes)


@dataclass
class ClientState:
    sock: socket.socket
    address: str
    molecule_id: int
    last_amp: Optional[np.ndarray] = None  # last source amplitude (3,)
    pending_send: bool = False
    initialized: bool = False
    alive: bool = True
    extras: dict = field(default_factory=dict)


def get_available_host_port() -> int:
    """Helper function to ask the OS for a local host and free TCP port.
    Returns (host, port) tuple.
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return "127.0.0.1", s.getsockname()[1]


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


class SocketHub:
    """
    A socket server that handles multiple driver connections with the FDTD engine.

    The major functionality of this class is to:

    - Accept many driver connections with the FDTD engine, including sending and receiving data.
    - Manage the lifecycle of driver connections, including initialization, polling, reconnection, and cleanup.
    """

    def __init__(
        self,
        host: Optional[str] = None,
        port: Optional[int] = 31415,
        unixsocket: Optional[str] = None,
        timeout: float = 60.0,
        latency: float = 0.01,
    ):
        """
        + **`host`** (str): Host address to bind the server. Default is None.
        + **`port`** (int): Port number to bind the server. Default is 31415.
        + **`unixsocket`** (str): Path to a Unix socket file. Default is None (no Unix socket).
        Using a Unix socket is recommended for **local connections** for the improved performance.
        If unixsocket is given a str, it will be created under /tmp/socketmxl_<str>. When using a
        Unix socket, the host and port arguments are ignored and should not be provided.
        + **`timeout`** (float): Timeout for socket operations in seconds. Default is 60.0 (s). This timeout controls
        when SocketHub will consider a driver connection to be dead and clean it up. If the computational time for
        each driver step is expected to be long, a larger timeout may be desirable.
        + **`latency`** (float): Latency for socket operations in seconds. Default is 0.01.
        For **local connections**, a lower latency, such as **1e-4**, may be desirable.
        """
        self.unixsocket_path = None
        if am_master():
            if unixsocket:
                self.serversock = socket.socket(socket.AF_UNIX)
                # mirror i-PI's /tmp/ipi_* default when given a name
                if not unixsocket.startswith("/"):
                    unixsocket = f"/tmp/socketmxl_{unixsocket}"
                self.unixsocket_path = unixsocket
                if os.path.exists(self.unixsocket_path):
                    probe = socket.socket(socket.AF_UNIX)
                    try:
                        probe.settimeout(0.25)
                        probe.connect(self.unixsocket_path)
                    except FileNotFoundError:
                        pass
                    except ConnectionRefusedError:
                        try:
                            os.unlink(self.unixsocket_path)
                        except FileNotFoundError:
                            pass
                    else:
                        probe.close()
                        raise RuntimeError(
                            f"Socket path {self.unixsocket_path} already in use"
                        )
                    finally:
                        try:
                            probe.close()
                        except Exception:
                            pass
                self.serversock.bind(unixsocket)
                self._where = unixsocket
            else:
                self.serversock = socket.socket(socket.AF_INET)
                self.serversock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                host = host or ""
                port = port or 31415
                self.serversock.bind((host, port))
                self._where = f"{host}:{port}"

            self.serversock.listen(4096)
            self.serversock.settimeout(0.25)

            self.timeout = float(timeout)
            self.latency = float(latency)

            # key: molecule_id or temp id
            self.clients: Dict[int, ClientState] = {}

            # peer -> molecule_id
            self.addrmap: Dict[str, int] = {}
            self._stop = False
            self._lock = threading.RLock()
            self._accept_th = threading.Thread(target=self._accept_loop, daemon=True)
            self._accept_th.start()

        # molecule_id -> ClientState (locked client)
        self.bound: Dict[int, ClientState] = {}

        # molecule ids we expect to serve
        self.expected: set[int] = set()

        # global pause when any driver is down
        self.paused = False

        # holds a frozen barrier until it successfully commits
        self._inflight = None

    def _accept_loop(self):
        """
        Accept loop thread: accept new connections and add to self.clients with temp id.
        """
        while not self._stop:
            try:
                csock, addr = self.serversock.accept()
            except socket.timeout:
                continue
            except OSError:
                break

            # NEW: trim latency and keep long-lived connections healthy
            try:
                csock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
                # Only for AF_INET; will raise on AF_UNIX
                csock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            except (OSError, AttributeError):
                pass  # AF_UNIX or platform without TCP_NODELAY

            peer = addr if isinstance(addr, str) else f"{addr[0]}:{addr[1]}"
            csock.settimeout(self.timeout)
            st = ClientState(sock=csock, address=peer, molecule_id=-1)
            with self._lock:
                # temp key: use id(csock) until INIT binds molecule_id
                self.clients[id(csock)] = st

    def _maybe_init_client(self, st: ClientState, init_payload: dict):
        """
        Handle NEEDINIT by sending INIT with JSON (contains molecule_id & unit scales).

        + **`st`** (ClientState): The client state to initialize.
        + **`init_payload`** (dict): The initialization payload to send.
        """
        pack_init(st.sock, init_payload)
        st.initialized = True

    def _dispatch_field(self, st: ClientState, efield_au: np.ndarray, meta: dict):
        """
        Send FIELDDATA (POSDATA) with the given E-field vector in a.u.

        + **`st`** (ClientState): The client state to send the field to.
        + **`efield_au`** (np.ndarray): The electric field vector in atomic units (a.u.).
        + **`meta`** (dict): Additional metadata to attach to this send.
        """
        try:
            send_msg(st.sock, FIELDDATA)
            I = np.eye(3, dtype=DT_FLOAT)
            send_array(st.sock, I.T, DT_FLOAT)
            send_array(st.sock, I.T, DT_FLOAT)
            send_int(st.sock, 1)
            vec = np.asarray(efield_au, dtype=DT_FLOAT).reshape(1, 3)
            send_array(st.sock, vec, DT_FLOAT)
            st.pending_send = True
            st.extras.update(meta or {})
        except (socket.timeout, SocketClosed, OSError):
            st.alive = False
            if st.molecule_id >= 0 and self.bound.get(st.molecule_id) is st:
                self._log(
                    f"DISCONNECTED (send): mol {st.molecule_id} from {st.address}"
                )
                self.bound[st.molecule_id] = None
            raise

    def _query_result(self, st: ClientState) -> Tuple[np.ndarray, bytes]:
        """
        Send GETSOURCE (GETFORCE) and read SOURCEREADY (FORCEREADY).

        + **`st`** (ClientState): The client state to query.

        Returns: (amplitude_vec3, extra_bytes)
        - amplitude_vec3: The source amplitude vector in the form [dmu_x/dt, dmu_y/dt, dmu_z/dt].
        - extra_bytes: Additional bytes sent by the driver.
        """
        try:
            send_msg(st.sock, GETSOURCE)
            msg = recv_msg(st.sock)
            if msg != SOURCEREADY:
                raise SocketClosed(f"Expected {SOURCEREADY!r}, got {msg!r}")
            e, forces, vir, extra = recv_getforce(st.sock)
            amp = np.array(forces[0], dtype=float)  # (3,)
            st.last_amp = amp
            st.pending_send = False
            return amp, extra
        except (socket.timeout, SocketClosed, OSError):
            st.alive = False
            if st.molecule_id >= 0 and self.bound.get(st.molecule_id) is st:
                self._log(
                    f"DISCONNECTED (recv): mol {st.molecule_id} from {st.address}"
                )
                self.bound[st.molecule_id] = None
            raise

    def _bind_client_locked(
        self, st: ClientState, molid: int, init_payload: dict, st_key
    ):
        """
        Bind client to molecule id if free; otherwise leave client unbound.

        + **`st`** (ClientState): The client state to bind.
        + **`molid`** (int): The molecule id to bind to.
        + **`init_payload`** (dict): The initialization payload to send if binding.
        + **`st_key`** (int): The temporary key of the client in self.clients.
        """
        if self.bound.get(molid) is None:
            self._maybe_init_client(st, init_payload)
            st.molecule_id = molid
            self.bound[molid] = st
            self.addrmap[st.address] = molid
            self.clients[molid] = st
            if st_key != molid:
                try:
                    del self.clients[st_key]
                except KeyError:
                    pass
            self._log(f"CONNECTED: mol {molid} <- {st.address}")
            # NEW: this molid is part of a frozen barrier -> force re-dispatch
            self._reset_inflight_for(molid)
            st.pending_send = False  # defensive: this is a fresh socket
            return True
        return False

    def _log(self, *a):
        """Log message helper"""
        print("[SocketHub]", *a)

    def _pause(self):
        """Pause the socket hub."""
        self.paused = True

    def _resume(self):
        """Resume the socket hub."""
        self.paused = False

    def _reset_inflight_for(self, molid: int):
        """
        If a barrier is frozen, force re-dispatch for this molid on reconnect.

        + **`molid`** (int): The molecule id to reset in the inflight barrier.
        """
        if self._inflight and (molid in self._inflight["wants"]):
            self._inflight["sent"][molid] = False
            self._inflight["ready"][molid] = False

    # -------------- public API --------------

    def register_molecule(self, molecule_id: int) -> None:
        """
        Reserve a slot for this molecule id (client may connect later).

        + **`molecule_id`** (int): The ID of the molecule.
        """
        with self._lock:
            # If already registered, raising a ValueError
            if molecule_id in self.expected:
                raise ValueError(f"Molecule ID {molecule_id} already registered!")
            # No explicit state needed yet; client binds on INIT.
            self.expected.add(int(molecule_id))
            self.bound.setdefault(int(molecule_id), None)

    def step_barrier(
        self, requests: Dict[int, dict], timeout: Optional[float] = None
    ) -> Dict[int, np.ndarray]:
        """
        Perform a barrier step: dispatch fields and collect source amplitudes. This method is the
        core of the SocketHub's functionality, coordinating the communication between the FDTD engine and
        multiple molecular drivers.

        + **`requests`** (dict): A dictionary mapping molecule IDs to their respective field requests.
          Each request should contain:
            - "efield_au": A tuple or list of three floats representing the electric field vector in atomic units (a.u.).
            - "meta": (Optional) A dictionary of additional metadata to attach to this send.
            - "init": (Optional) A dictionary containing initialization parameters for the molecule.
        + **`timeout`** (float, optional): Maximum time to wait for the barrier step to complete.
        If None, uses the default timeout set during initialization.
        """
        if self.paused:
            return {}

        deadline = time.time() + (timeout or self.timeout)
        results: Dict[int, dict] = {}

        # If a barrier is already in flight, ignore new 'requests' and reuse the frozen one.
        if self._inflight is None:
            wants = set(requests.keys())
            self._inflight = {
                "wants": wants,
                "efields": {
                    mid: np.asarray(requests[mid]["efield_au"], dtype=DT_FLOAT).copy()
                    for mid in wants
                },
                "meta": {mid: requests[mid].get("meta", {}) for mid in wants},
                "sent": {mid: False for mid in wants},
                "ready": {mid: False for mid in wants},
            }
        else:
            # Reuse the frozen barrier even if the caller passed different fields
            wants = set(self._inflight["wants"])

        ready = self._inflight["ready"]

        # --- hard gate: do not dispatch fields until everyone is bound ---
        ids = set(int(k) for k in requests.keys())
        with self._lock:
            if not self.all_bound(ids, require_init=True):
                # Try to progress INIT quickly, but DO NOT send FIELDDATA yet
                for st_key, st in list(self.clients.items()):
                    if not st or not st.alive:
                        continue
                    try:
                        send_msg(st.sock, STATUS)
                        reply = recv_msg(st.sock)
                    except (socket.timeout, SocketClosed, OSError):
                        st.alive = False
                        if st.molecule_id >= 0 and self.bound.get(st.molecule_id) is st:
                            self._log(
                                f"DISCONNECTED: mol {st.molecule_id} from {st.address}"
                            )
                            self.bound[st.molecule_id] = None
                            self._pause()
                            # NEW: make the frozen barrier re-send the old field to this molid
                            self._reset_inflight_for(st.molecule_id)
                        continue
                    if reply == NEEDINIT:
                        for mid in ids:
                            if self.bound.get(mid) is None:
                                init_payload = requests.get(mid, {}).get(
                                    "init", {"molecule_id": mid}
                                )
                                self._bind_client_locked(
                                    st, int(mid), init_payload, st_key
                                )
                                break
                return {}  # nothing dispatched; drivers remain idle

        # --- normal step_barrier continues below ---
        aborted = False
        with self._lock:
            # 1. poll and init/dispatch
            for st_key, st in list(self.clients.items()):
                try:
                    send_msg(st.sock, STATUS)
                    reply = recv_msg(st.sock)
                except (socket.timeout, SocketClosed, OSError):
                    st.alive = False
                    if st.molecule_id >= 0 and self.bound.get(st.molecule_id) is st:
                        self._log(
                            f"DISCONNECTED: mol {st.molecule_id} from {st.address}"
                        )
                        self.bound[st.molecule_id] = None
                        self._pause()
                        # NEW: make the frozen barrier re-send the old field to this molid
                        self._reset_inflight_for(st.molecule_id)
                    aborted = True
                    continue
                if reply == NEEDINIT:
                    # Bind this client to the first UNBOUND molecule id present in requests
                    # Preserve the requests' order (Python dicts keep insertion order).
                    chosen = None
                    for mid in requests.keys():
                        if self.bound.get(mid) is None:
                            chosen = mid
                            break
                    # Fallback: allow a rebind if this client was previously bound and crashed/reconnected
                    if (
                        chosen is None
                        and st.molecule_id >= 0
                        and self.bound.get(st.molecule_id) is None
                    ):
                        chosen = st.molecule_id
                    if chosen is None:
                        # Nothing to serve right now; keep the client idle
                        continue
                    init_payload = requests.get(chosen, {}).get(
                        "init", {"molecule_id": chosen}
                    )
                    self._bind_client_locked(st, int(chosen), init_payload, st_key)

                elif reply == READY:
                    molid = st.molecule_id
                    if (
                        molid in self._inflight["wants"]
                        and not st.pending_send
                        and not self._inflight["sent"][molid]
                    ):
                        evec = self._inflight["efields"][molid]
                        meta = self._inflight["meta"][molid]
                        try:
                            self._dispatch_field(st, evec, meta)
                            self._inflight["sent"][molid] = True
                        except Exception:
                            aborted = True
                            break

                elif reply == HAVEDATA:
                    if st.molecule_id in ready:
                        ready[st.molecule_id] = True
                    continue

            # 2. second pass: wait for all pending molecules to finish (barrier)
            wants = set(self._inflight["wants"])

        if aborted:
            # abort this barrier; caller will enter the pause path
            return {}

        while time.time() < deadline and (wants - set(results.keys())):
            time.sleep(self.latency)

            with self._lock:
                # Always sweep over ALL clients, not just molecule-ids in 'wants'
                for st_key, st in list(self.clients.items()):
                    if not st or not st.alive:
                        continue
                    try:
                        send_msg(st.sock, STATUS)
                        reply = recv_msg(st.sock)
                    except (socket.timeout, SocketClosed, OSError):
                        st.alive = False
                        if st.molecule_id >= 0 and self.bound.get(st.molecule_id) is st:
                            self._log(
                                f"DISCONNECTED: mol {st.molecule_id} from {st.address}"
                            )
                            self.bound[st.molecule_id] = None
                            self._pause()
                            # NEW: make the frozen barrier re-send the old field to this molid
                            self._reset_inflight_for(st.molecule_id)
                        continue

                    if reply == NEEDINIT:
                        # Bind this client to the first UNBOUND molecule id present in requests
                        # Preserve the requests' order (Python dicts keep insertion order).
                        chosen = None
                        for mid in requests.keys():
                            if self.bound.get(mid) is None:
                                chosen = mid
                                break
                        # Fallback: allow a rebind if this client was previously bound and crashed/reconnected
                        if (
                            chosen is None
                            and st.molecule_id >= 0
                            and self.bound.get(st.molecule_id) is None
                        ):
                            chosen = st.molecule_id
                        if chosen is None:
                            # Nothing to serve right now; keep the client idle
                            continue
                        init_payload = requests.get(chosen, {}).get(
                            "init", {"molecule_id": chosen}
                        )
                        self._bind_client_locked(st, int(chosen), init_payload, st_key)
                        continue  # next client

                    if reply == READY:
                        molid = st.molecule_id
                        if (
                            molid in self._inflight["wants"]
                            and not st.pending_send
                            and not self._inflight["sent"][molid]
                        ):
                            evec = self._inflight["efields"][molid]
                            meta = self._inflight["meta"][molid]
                            self._dispatch_field(st, evec, meta)
                        continue

                    if reply == HAVEDATA:
                        if st.molecule_id in ready:
                            ready[st.molecule_id] = True
                        continue

                    # Some clients may erroneously send STATUS; ignore gracefully
                    if reply == STATUS:
                        continue

                # exit condition: filled everything we wanted
                if all(ready.get(mid, False) for mid in wants):
                    break

        # Abort if not everyone is ready (e.g., disconnect); keep the frozen barrier.
        if not all(
            self._inflight["ready"].get(mid, False) for mid in self._inflight["wants"]
        ):
            return {}

        # Phase C: commit all together (send GETSOURCE to everyone now)
        for mid in self._inflight["wants"]:
            st = self.clients.get(mid)
            if not st or not st.alive:
                return {}  # keep barrier for retry
            try:
                amp, extra = self._query_result(st)
                results[mid] = {"amp": amp, "extra": extra}
            except (socket.timeout, SocketClosed, OSError):
                return {}  # keep barrier for retry

        # SUCCESS — clear the frozen barrier
        self._inflight = None

        # Fallback: any missing results -> use last known amplitude (if any),
        # We turn off this fallback for now to avoid silent errors.
        if False:
            with self._lock:
                for mid in wants - set(results.keys()):
                    st = self.clients.get(mid)
                    if st and st.last_amp is not None:
                        results[mid] = {"amp": st.last_amp.copy(), "extra": b""}
        return results

    def all_bound(self, molecule_ids, require_init=True):
        """
        Check if all given molecule_ids are bound (and optionally initialized).

        + **`molecule_ids`** (iterable): An iterable of molecule IDs to check.
        + **`require_init`** (bool): If True, also require that the clients are initialized. Default is True.
        """
        with self._lock:
            for mid in molecule_ids:
                st = self.bound.get(int(mid))
                if st is None or not st.alive:
                    return False
                if require_init and not st.initialized:
                    return False
            return True

    def wait_until_bound(self, init_payloads: dict, require_init=True, timeout=None):
        """
        BLOCK until all molecule_ids in init_payloads are bound.

        + **`init_payloads`** (dict): A dictionary mapping molecule IDs to their respective initialization payloads.
          Each payload should be a dictionary containing initialization parameters for the molecule.
        + **`require_init`** (bool): If True, also require that the clients are initialized. Default is True.
        + **`timeout`** (float, optional): Maximum time to wait for all clients to bind. If None, uses the default timeout set during initialization.
        """
        wanted = {int(k) for k in init_payloads.keys()}
        deadline = time.time() + (timeout if timeout is not None else self.timeout)

        while True:
            if self.all_bound(wanted, require_init=require_init):
                self._resume()
                return True

            # Progress the INIT handshakes without dispatching any field steps
            with self._lock:
                # Only touch clients for molecule-ids that are still unbound
                pending_ids = {mid for mid in wanted if self.bound.get(mid) is None}
                if not pending_ids:
                    # (shouldn't happen because of the all_bound check, but be safe)
                    time.sleep(self.latency)
                    continue

                for st_key, st in list(self.clients.items()):
                    if not st or not st.alive:
                        continue
                    # Only ping clients that are currently unbound (either brand-new or reconnecting)
                    if st.molecule_id >= 0 and st.molecule_id not in pending_ids:
                        continue
                    try:
                        send_msg(st.sock, STATUS)
                        reply = recv_msg(st.sock)
                    except (socket.timeout, SocketClosed, OSError):
                        st.alive = False
                        # free a binding if this was a re-connect case
                        if st.molecule_id >= 0 and self.bound.get(st.molecule_id) is st:
                            self._log(
                                f"DISCONNECTED: mol {st.molecule_id} from {st.address}"
                            )
                            self.bound[st.molecule_id] = None
                            self._pause()
                            # NEW: make the frozen barrier re-send the old field to this molid
                            self._reset_inflight_for(st.molecule_id)
                        continue

                    if reply == NEEDINIT:
                        # choose an unclaimed id from 'wanted'
                        chosen = None
                        for mid in pending_ids:
                            if self.bound.get(mid) is None:
                                chosen = mid
                                break
                        if chosen is not None:
                            self._bind_client_locked(
                                st, int(chosen), init_payloads[int(chosen)], st_key
                            )

            if timeout is not None and time.time() > deadline:
                return False
            time.sleep(self.latency)

    def graceful_shutdown(self, reason: Optional[str] = None, wait: float = 2.0):
        """
        Politely ask all connected drivers to exit, and wait briefly for BYE.

        + **`reason`** (str): Optional reason for shutdown to log. Default is None.
        + **`wait`** (float): Time in seconds to wait for drivers to respond. Default is 2.0 (s).
        """
        with self._lock:
            for st in list(self.clients.values()):
                if not st or not st.alive:
                    continue
                try:
                    send_msg(st.sock, STOP)
                except Exception:
                    st.alive = False
                    if st.molecule_id >= 0 and self.bound.get(st.molecule_id) is st:
                        self._log(
                            f"DISCONNECTED: mol {st.molecule_id} from {st.address}"
                        )
                        self.bound[st.molecule_id] = None
                        self._pause()

        deadline = time.time() + float(wait)
        while time.time() < deadline:
            time.sleep(self.latency)
            with self._lock:
                for st in list(self.clients.values()):
                    if not st or not st.alive:
                        continue
                    try:
                        # Make reads snappy during shutdown
                        st.sock.settimeout(self.latency)
                        msg = recv_msg(st.sock)
                        if msg == BYE:
                            # Clean close on our side
                            st.alive = False
                            if (
                                st.molecule_id >= 0
                                and self.bound.get(st.molecule_id) is st
                            ):
                                self._log(
                                    f"DISCONNECTED: mol {st.molecule_id} from {st.address}"
                                )
                                self.bound[st.molecule_id] = None
                            try:
                                st.sock.shutdown(socket.SHUT_RDWR)
                            except Exception:
                                pass
                            try:
                                st.sock.close()
                            except Exception:
                                pass
                    except (socket.timeout, SocketClosed, OSError):
                        # Either no message yet or peer closed already; keep sweeping
                        continue

    def stop(self):
        """
        Stop accept loop, tell clients to exit, then close everything.
        """
        # First, stop accepting new connections
        self._stop = True
        try:
            self.serversock.close()
        except Exception:
            pass

        # Then, gracefully end existing sessions
        try:
            self.graceful_shutdown(wait=max(2.0, 10 * self.latency))
        finally:
            with self._lock:
                for st in list(self.clients.values()):
                    try:
                        st.sock.close()
                    except Exception:
                        pass

        # if unix socket, remove the path
        if self.unixsocket_path and os.path.exists(self.unixsocket_path):
            os.unlink(self.unixsocket_path)
            print(f"[SocketHub] Unlinked unix socket path {self.unixsocket_path}")
