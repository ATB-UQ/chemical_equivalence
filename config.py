from os.path import exists, abspath, join, dirname
import yaml

DOUBLE_BOND_LENGTH_CUTOFF = {
    frozenset(['C', 'C']): 0.1380, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "C", 1.5)]
    frozenset(['C', 'N']): 0.1337, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "N", 1.5)]
    frozenset(['N', 'N']): 0.1250, #nm, Source: http://www.chemikinternational.com/wp-content/uploads/2014/04/13.pdf
}

THIS_DIR = dirname(abspath(__file__))
YAML_CONFIG_FILE = join(THIS_DIR, "config.yml")

if not exists(YAML_CONFIG_FILE):
    raise Exception("config.yml file not found, modify config.yml.example file change it's name to config.yml")

with open(YAML_CONFIG_FILE) as fp:
    yaml_config = yaml.load(fp)

NAUTY_EXECUTABLE = yaml_config["NAUTY_EXECUTABLE"] \
    if yaml_config["NAUTY_EXECUTABLE"].startswith("/") \
    else join(THIS_DIR, yaml_config["NAUTY_EXECUTABLE"])

assert NAUTY_EXECUTABLE[0] == '/', 'Dreadnaut executable was not an absolute path: "{0}"'.format(NAUTY_EXECUTABLE)

assert exists(NAUTY_EXECUTABLE), 'Could not find dreadnaut executable at: "{0}". Did you install nauty (http://users.cecs.anu.edu.au/~bdm/nauty/) ?'.format(NAUTY_EXECUTABLE)
