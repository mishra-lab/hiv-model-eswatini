import os,re
docs = ['intro','methods','results','discuss']
igs = [
  '(?<!\\\\)%.*',
  '\\\\par',
  '\\\\item',
# '\\s*?\\\\cite\{.*?\}',
  '\\\\input\{.*?\}',
  '\\\\label\{.*?\}',
  '\\\\(begin|end)\{(enumerate|itemize)\}',
  re.compile('\\\\begin\{table\}(.*?)\\\\end\{table\}',re.DOTALL),
  re.compile('\\\\begin\{figure\}(.*?)\\\\end\{figure\}',re.DOTALL),
]
reps = [
  ('~', ' '),
  ('``','"'),
  ('\'\'','"'),
  ('\\$(.*?)\\$', '\\1'),
  ('\\\\times', 'x'),
  ('\\\\emph\{(.*?)\}', '\\1'),
  ('\\\\text..\{(.*?)\}', '\\1'),
  ('\\\\(?:sub)*section\{(.*?)\}', '\n\\1\n'),
  ('\\\\ref\{.*?\}', 'REF'),
  ('\\\\x\{.*?\}', 'X'),
]

body = ''

for doc in docs:
  with open(doc+'.tex','r') as f:
    body += '\n\n'+doc.upper()+'\n\n'+f.read()

for k,v in reps:
  body = re.sub(k,v,body)

for i in igs:
  body = re.sub(i,'',body)

for n in range(3):
  body = body.replace('\n\n\n','\n\n')

with open('wc/body.tmp','w') as f:
  f.write(body)

os.system('wc -w wc/body.tmp | cut -d " " -f1 && rm wc/body.tmp')
os.system('wc -w abstract.tex | cut -d " " -f1')

