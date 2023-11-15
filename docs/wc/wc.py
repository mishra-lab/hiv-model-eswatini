import os,sys,re
# setup
fins = sys.argv[1:]
igs = [
  '(?<!\\\\)%.*',
  '\\\\(par|item)\\s',
  # '\\s*?\\\\cite\{.*?\}',
  # '\\\\input\{.*?\}',
  '\\\\label\{.*?\}',
  '\\\\(begin|end)\{(enumerate|itemize)\}',
]
envs = [
  re.compile('\\s*?\\\\begin\{'+env+'\}(.*?)\\\\end\{'+env+'\}',re.DOTALL)
  for env in ['table','figure','equation','alignat']
]
reps = [
  ('~', ' '),
  ('(``|\'\')', '"'),
  ('\\$(.*?)\\$', '\\1'),
  ('\\\\times', 'x'),
  ('\\\\emph\{(.*?)\}', '\\1'),
  ('\\\\text..\{(.*?)\}', '\\1'),
  ('\\\\(?:sub)*section\{(.*?)\}', '\n\\1\n'),
  ('\\\\paragraph\{(.*?)\}', '\n\\1\n'),
  ('\\\\ref\{(.*?)\}', '\\1'),
  ('\\\\cite\{(.*?)\}', '[\\1]'),
  ('\\\\href\{.*?\}\{(.*?)\}', '\\1'),
]
# load
body = ''
for fin in fins:
  with open(fin+'.tex','r') as f:
    body += '\n\n<<<'+fin.upper()+'>>>\n\n'+f.read()
# clean
for k,v in reps:   body = re.sub(k,v,body)
for i in igs+envs: body = re.sub(i,'',body)
for n in range(9): body = body.replace('\\s*\n','\n').replace('\n\n\n','\n\n')
# output
with open('body.tmp','w') as f:
  f.write(body)
os.system('echo -n body:\  && wc -w body.tmp | cut -d " " -f1 && rm body.tmp')
os.system('echo -n abs:\ \ && wc -w abstract.tex | cut -d " " -f1')
