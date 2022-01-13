# This Python file uses the following encoding: utf-8
import strictyaml as syaml
class Microscope:
    yaml = None
    def __init__(self):
        pass
    def __init__(self, path):
        yaml = syaml.load(path)

    def test_connection(self):
        return False
