"""Base objects for ps"""

from functools import partial
from typing import Iterable, Callable, Optional, Union

# Note: because typing.Mapping made Commands have signature (*args, **kwargs):
from collections.abc import Mapping
from ps.util import (
    run,
    str_if_bytes,
    local_identifier_command_dict,
    IdentifierCommandDict,
    IdentifiedCommands,
)

from dataclasses import dataclass


def mk_raw_command_func(command: str):
    def run_command(command_args_str: str = ''):
        return run_command(command_args_str + ' ' + command)

    return run_command


def join_if_not_string(iterable: Iterable) -> str:
    if not isinstance(iterable, str):
        iterable = ' '.join(iterable)
    return iterable


# TODO: Figure out the standard CLI language: command, instruction, argument, parameter


def first_valid_result(
    funcs, *args, _is_valid=bool, _ignore_exceptions=(Exception,), **kwargs
):
    for func in funcs:
        try:
            result = func(*args, **kwargs)
            if _is_valid(result):
                return result
        except _ignore_exceptions:
            pass


def man_1_page_str(command):
    return str_if_bytes(run(f'man 1 {command}'))


def dash_dash_help_str(command):
    return str_if_bytes(run(f'{command} --help'))


get_doc_options = {
    'man_1_page_str': man_1_page_str,
    'dash_dash_help_str': dash_dash_help_str,
}


def find_doc(command, *, doc_finders=(man_1_page_str, dash_dash_help_str)):
    return first_valid_result(doc_finders, command)


_dflt_run = run  # to use when context overwrites "run" name


class Command:
    """
    A ``Command`` runs a specific shell script for you in a specific manner.
    The ``run`` function is the general function to do that, and we saw
    that you can curry `run` to specify what and how to run.
    ``Command`` just wraps such a curried ``run`` function
    (or any compliant run function you provide),
    and specifies what executable (the ``command`` argument) to actually run.

    So not much over a curried ``run``.

    But what it does do as well is set up the ability to do other things
    that may be specific to the executable you're running, such as
    giving your (callable) command instance a signature, some docs, or
    a help method.

    >>> import os
    >>> pwd = Command('pwd')
    >>> os.path.isdir(pwd())
    True
    >>> assert pwd.__doc__  # docs exist (and are non-empty)!
    >>> # To print the docs:
    >>> pwd.help()  # doctest: +SKIP
    PWD(1)                    BSD General Commands Manual                   PWD(1)

    NAME
         pwd -- return working directory name

    SYNOPSIS
         pwd [-L | -P]
    ...


    """

    def __init__(
        self,
        command: Union[str],
        run: Optional[Callable] = None,
        get_doc: Callable[[str], str] = find_doc,
        **run_kwargs,
    ):
        self.command = command
        run = run or _dflt_run
        if run_kwargs:
            run = partial(run, **run_kwargs)
        self.run = run
        self.get_doc = get_doc

    def __call__(self, args: Iterable = ''):
        return self.raw_call(args)

    def raw_call(self, args: Iterable = ''):
        return self.run(self.instruction_str(args))

    def instruction_str(self, args: Iterable):
        args_str = join_if_not_string(args)
        return f'{self.command} {args_str}'

    @property
    def __doc__(self):
        return self.help_str()

    # TODO: Include 'intelligence' to find the appropriate help string
    def help_str(self, get_doc=None):
        get_doc = get_doc or self.get_doc
        if isinstance(get_doc, str):
            if get_doc not in get_doc_options:
                raise ValueError(
                    f'get_doc not found. Specify a function or one of {get_doc_options}'
                )
            get_doc = get_doc_options[get_doc]
        return get_doc(self.command) or ''

    def help(self, get_doc=None):
        return print(self.help_str(get_doc))


def _ensure_identifier_keyed_dict(
    commands: IdentifiedCommands,
) -> IdentifierCommandDict:
    if callable(commands):
        commands = commands()
    if isinstance(commands, Mapping):
        d = dict(commands)
    else:
        d = {x: x for x in commands}
    non_identifier_keys = list(filter(lambda x: not x.isidentifier(), d))
    if non_identifier_keys:
        raise ValueError(f'These were not identifiers: {non_identifier_keys}')
    return d


class Commands(Mapping):
    r"""
    A collection of commands.

    The general usage is that you can specify a mapping between valid python identifiers
    (alphanumerical strings (and underscores) that don't start with a number) and
    functions. If instead of functions you specify a string, a ``factory`` comes
    in play to make a function based on your string.
    By default, it will consider it as a console command and give you a function that
    runs it.

    >>> import os
    >>>
    >>> c = Commands({
    ...     'current_dir': 'pwd',
    ...     'sys_listdir': 'ls -l',
    ...     'listdir': os.listdir,
    ...     'echo': Command('echo', egress=lambda x: print(x.decode().strip())),
    ... })
    >>>
    >>> list(c)
    ['current_dir', 'sys_listdir', 'listdir', 'echo']
    >>> current_folder = c.current_dir()
    >>> os.path.isdir(current_folder)
    True
    >>> b = c.sys_listdir()
    >>> b[:40]  # doctest: +SKIP
    b'total 56\n-rw-r--r--@ 1 Thor.Whalen  staf'
    >>> a_list_of_filenames = c.listdir()
    >>> isinstance(a_list_of_filenames, list)
    True
    >>> c.echo('hello world')
    hello world

    If you don't specify any commands, it will gather all executable names it can find in
    your local system (according to your ``PATH`` environment variable),
    map those to valid python identifiers if needed, and use that.

    Important: Note that finding executable in the ``PATH`` doesn't mean that it will
    work, or is safe -- so use with care!

    >>> c = Commands()
    >>> assert len(c) > 0
    >>> 'ls' in c
    True

    You can access the 'ls' command as a key or an attribute

    >>> assert c['ls'] == c.ls

    You can print the ``.help()`` (docs) of any command, or just get the help string:

    >>> man_page = c.ls.help_str()

    Let's see if these docs have a few things we expect it to have for ``ls``:

    >>> assert 'BSD General Commands Manual' in man_page  # doctest: +SKIP
    >>> assert 'list directory contents' in man_page

    Let's see what the output of ``ls`` gives us:

    >>> output = c.ls('-la').decode()  # execute "ls -la"
    >>> assert 'total' in output  # "ls -l" output usually starts with "total"
    >>> assert '..' in output  # the "-a" flag includes '..' as a file

    Note that we needed to decode the output here.
    That's because by default the output of a command will be captured in bytes.
    If you want to apply a decoder to (attempt to) convert all outputs into strings,
    you can specify a ``factory`` that will do this for you automatically.

    The default ``factory`` is ``Command``, which  has a ``run``
    argument that defines how an instruction should be run. 
    The default of ``run`` is ``run_command``, which conveniently has an ``egress`` 
    argument where you can specify a function to call on the output. 
    
    So one solution to define a ``Commands`` instance that will (attempt to) output 
    strings systematically is to do this:

    >>> from functools import partial
    >>> from ps.util import run
    >>> run_and_cast_to_str = partial(run, egress=bytes.decode)
    >>> factory = partial(Command, run=run_and_cast_to_str)
    >>> cc = Commands(factory=factory)

    So now we have:

    >>> output = cc.ls('-la')
    >>> isinstance(output, str)  # it's already decoded for us!
    True

    """

    def __init__(
        self,
        commands: IdentifiedCommands = local_identifier_command_dict,
        factory: Callable[[str], Callable] = Command,
    ):
        self._commands = _ensure_identifier_keyed_dict(commands)
        for name, command in self._commands.items():
            if callable(command):
                # if command value is already a callable we just take it as is:
                runner = command
            else:
                runner = factory(command)
            setattr(self, name, runner)
        self._factory = factory

    def __getitem__(self, command):
        try:
            return getattr(self, command)
        except AttributeError:
            raise KeyError(f"This command was not found: '{command}'")

    def __iter__(self):
        yield from self._commands

    def __len__(self):
        return len(self._commands)

    def __contains__(self, command):
        return command in self._commands
