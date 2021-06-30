from growth import *
from pyflowchart import Flowchart

with open('Main.py') as f:
        code = f.read()

fc = Flowchart.from_code(code, field='', inner=True, simplify=True)
print(fc.flowchart())




