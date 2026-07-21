# Standard Library Imports
from datetime import datetime
from typing import Optional

# Third-Party Imports

# Local Library Imports
from consts.dtypes import NP_INT

def print_msg(msg: str, l_rank: NP_INT, print_rank: Optional[NP_INT] = None):
    current_time: str = datetime.now().strftime("%H:%M:%S")
    
    out_msg: str
    if print_rank is None:
        out_msg = "[{}]: [{}]: {}".format(current_time, l_rank, msg)
        print(out_msg, flush = True)
    elif l_rank == print_rank:
        out_msg = "[{}]: {}".format(current_time, msg)
        print(out_msg, flush = True)