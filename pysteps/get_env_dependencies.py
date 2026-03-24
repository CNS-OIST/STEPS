import os

__all__ = ["dynamic_metadata"]


def __dir__():
    return __all__


KEYS = {"env-prefix"}

DEFAULT_EXTRA_PREFIX = "default-"


def dynamic_metadata(field, settings, project):
    """Return the list of package dependencies
    This function reads dependencies from the project.optional-dependencies table
    All extras that start with DEFAULT_EXTRA_PREFIX are added to the list, unless a
    specific environment variable is present and set to true.
    """
    if settings.keys() > KEYS:
        msg = f"Only {KEYS} settings allowed by this plugin"
        raise RuntimeError(msg)

    if "env-prefix" not in settings:
        msg = "Must contain the 'env-prefix' setting"
        raise RuntimeError(msg)

    all_deps = []
    for extra, deps in project["optional-dependencies"].items():
        if extra.startswith(DEFAULT_EXTRA_PREFIX):
            extra_name = extra[len(DEFAULT_EXTRA_PREFIX) :]
            env_var = settings["env-prefix"] + extra_name.upper()
            # Only add the dependencies if there is no env variable preventing it
            if os.environ.get(env_var, "false").lower() != "true":
                all_deps += deps

    return all_deps
