"""

"""

from ps.base import Command, Commands

_locals = locals()
for _command, _command_func in Commands().items():
    _locals[_command] = _command_func

del _locals
del _command
del _command_func
