verb = 3

def log(lvl,msg=''):
  if lvl < 0: print(''); return msg # flush log
  pre = ['-'*80+'\n','',' > ',''][lvl]
  end = ['\n','\n','\n',''][lvl]
  if lvl <= verb:
    print(pre+str(msg),end=end,flush=True)
