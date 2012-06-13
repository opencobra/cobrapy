# attempt to import all working solvers in this directory
solver_list = {}
from os import listdir as _listdir
from os import path as _path
for i in _listdir(_path.split(_path.abspath(__file__))[0]):
    if i.startswith("_") or i.startswith("."):
        continue
    if not i.endswith(".py"):
        continue
    try:
        exec("import %s" % i.strip(".py"))
        solver_list[i.strip(".py")] = eval(i.strip(".py"))
    except Exception, e:
        print Exception, e
        pass
del _path
del _listdir
del i
