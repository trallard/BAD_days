from os import remove, close
import re, glob, os
from shutil import move
from tempfile import mkstemp

## extra needed to add '' in the class or it fails to render
rdict = {'class=': 'class="', '>':'">', 'scope=row':'scope="row','scope=col':'scope="col'}
s_repl = {'class=', 'scope=row', 'scope=col'}
robj = re.compile('|'.join(rdict.keys()))


# replace the expressions that are not correcly parsed
def replace(file_path):
    #Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path,'w') as new_file:
        with open(file_path) as source_file:
            for line in source_file:
                if any(s in line for s in s_repl):
                    line_new = robj.sub(lambda m: rdict[m.group(0)], line)
                    new_file.write(line_new)
                    print(line)
                else:
                    new_file.write(line)
    new_file.close()
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)


# getting the notebooks in the appropriate file
def locate(nb_dir):
    abs_path = os.path.join('.', nb_dir, '*.md')
    notebooks = glob.glob(abs_path)
    return(notebooks)


notebooks = locate('notebooks')
for nb in notebooks: replace(nb)
