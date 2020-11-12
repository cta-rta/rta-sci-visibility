# Visibility tools for rta-sci

## Installation

Create a virtual environment with all required dependences: 
```bash
conda env create --name <envname> --file=env_visibility.yaml
```

### Execution

Modify the configuration file:
```yaml
path:
    filename: <path/to/input/fits/file>
...
```

Run the available example script:
```bash
python visibility_test.py config_visibility.yaml
```
