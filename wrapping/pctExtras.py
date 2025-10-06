import itk
from itk import PCT as pct
import importlib
import shlex


def make_application_func(app_name):
    """
    Factory: returns a function that runs the PCT application `app_name`
    with either Python-style kwargs or a single CLI-style string.
    """

    def app_func(*args, **kwargs):
        app_module = importlib.import_module(f"itk.{app_name}")
        parser = app_module.build_parser()
        # Ensure help/usage shows the logical app name in Python contexts
        parser.prog = app_name

        if kwargs:
            if hasattr(parser, "parse_kwargs"):
                args_ns = parser.parse_kwargs(func_name=app_name, **kwargs)
            else:
                raise TypeError(f"Parser for {app_name} has no parse_kwargs method.")
        elif args and len(args) == 1 and isinstance(args[0], str):
            # Treat single string argument as CLI
            argv = shlex.split(args[0])
            args_ns = parser.parse_args(argv)
        else:
            args_ns = parser.parse_args()

        return app_module.process(args_ns)

    # Metadata for help()
    _parser = importlib.import_module(f"itk.{app_name}").build_parser()
    _parser.prog = app_name
    _parser.apply_signature(app_func)
    app_func.__name__ = app_name
    app_func.__module__ = "itk.PCT"
    # Python-only help: version, description + examples + options header
    description = (_parser.description or "").rstrip()
    examples = _parser.build_usage_examples(app_name)
    options = _parser.format_help()
    idx = options.lower().find("options:")
    opt_text = options[idx:].strip()
    parts = [pct.version(), description, examples, opt_text]
    app_func.__doc__ = "\n\n".join(parts)

    return app_func
