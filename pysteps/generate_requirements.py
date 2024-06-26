###___license_placeholder___###

import os
import sys
try:
    import tomllib
except ImportError:
    # TODO: Drop the dependency to toml when we require python >=3.11
    import toml as tomllib

if __name__ == '__main__':
    # This script parses pyproject.toml and generates a requirement file.
    #
    # It takes the following command line arguments:
    #   - A path to the pyproject.toml file to parse
    #   - A path to the requirement file to write
    #   - A comma-separated list of "extras" (optional)
    # See https://packaging.python.org/en/latest/guides/writing-pyproject-toml/
    # for more information about extras in pyproject.toml files.
    #
    # If something fails, exit with exit code 1.
    if not 3 <= len(sys.argv) <= 4:
        print(
            f'Error: expected 2 or 3 command line arguments, got {len(sys.argv)-1}: {sys.argv[1:]}', file=sys.stderr)
        sys.exit(1)
    toml_path = sys.argv[1]
    requirement_path = sys.argv[2]
    if len(sys.argv) > 3:
        opt_list = list(map(str.strip, sys.argv[3].split(',')))
    else:
        opt_list = []
    if not os.path.isfile(toml_path):
        print(
            f'Error: could not find file TOML file: {toml_path}', file=sys.stderr)
        sys.exit(1)
    with open(toml_path, 'r') as f:
        proj = tomllib.loads(f.read())
        project = proj.get('project', None)
        if project is None:
            print('No project section in the TOML file.', file=sys.stderr)
            sys.exit(1)
        deps = set(project.get('dependencies', []))
        opt_deps = project.get('optional-dependencies', {})
        for opt in opt_list:
            deps.update(opt_deps.get(opt, []))
    with open(requirement_path, 'w') as f:
        f.write('\n'.join(sorted(deps)))
