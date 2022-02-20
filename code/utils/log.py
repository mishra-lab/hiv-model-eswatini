verb = 3

def log(lvl,msg=''):
  pre = ['-'*80+'\n','',' > ',''][lvl]
  end = ['\n','\n','\n',''][lvl]
  if lvl <= verb:
    print(pre+str(msg),end=end,flush=True)
  