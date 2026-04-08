from casampi.MPICommandClient import MPICommandClient, MPIEnvironment
from traceback import print_exc
MPICLIENT = None


def get_mpi_client():
    """Returns the existing client or starts a new one if needed."""
    global MPICLIENT
    if MPIEnvironment.is_mpi_enabled:
        if MPICLIENT is None:
            start_mpi()
    return MPICLIENT

def start_mpi():
    global MPICLIENT
    try:
        MPICLIENT = MPICommandClient()
        MPICLIENT.start_services()
        MPICLIENT.set_log_mode('redirect')
        MPICLIENT.set_log_level('INFO4')
    except Exception as e:
        print_exc()
        raise
    return MPICLIENT
    