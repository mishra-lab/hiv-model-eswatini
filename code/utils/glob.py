_glob = {}

def globs(get=None,clear=None,iter=False,default=None,**kwds):
  # tool for getting, clearing, and setting global variables (one dictionary)
  # get:     key to get
  # clear:   key to clear
  # **kwds:  "key=value" variable assignments
  # iter:    for kwds: is value iterable & should we append to it, vs overwrite?
  # default: for get & clear: default get value if ket not found
  global _glob
  if get:
    return _glob.get(get,default)
  if clear:
    return _glob.pop(clear,default)
  for key,item in kwds.items():
    if iter:
      if key in _glob:
        _glob[key].append(item)
      else:
        _glob[key] = [item]
    else:
      _glob[key] = item
  return None