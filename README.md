# Visibility tools for rta-sci

## Installation

Create a virtual environment with all required dependencies: 
```bash
conda env create --name <envname> --file=env_visibility.yaml
```

### Example

Modify the configuration file (es: config_visibility.yaml):
```yaml
path:
    filename: <path/to/input/fits/file>
...
```

Run the available example script:
```bash
python visibility_test.py config_visibility.yaml
```
