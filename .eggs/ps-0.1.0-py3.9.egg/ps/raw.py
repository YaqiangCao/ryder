"""
Access to the raw commands of the system.

That is, you will be specifying the arguments of the commands through a single string.

"""

from ps.util import local_identifier_command_dict as _get_command_for_id
from ps.base import Command as _Command

for _id, _command in _get_command_for_id().items():
    globals()[_id] = _Command(_command)
