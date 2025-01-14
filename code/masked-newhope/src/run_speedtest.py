from os import system, popen
import sys
PATH = "./bench_res.txt"
open(PATH, 'w').close() # clear file
cur = ""
res = ""

VERBOSE_COMPILE = True
REDIRECT = ""
MAX_ORDER = 9


if not VERBOSE_COMPILE:
	REDIRECT=">/dev/null"
with open(PATH,'a') as f:
  for rng in range(1,2):
    for i in range(1, MAX_ORDER+1):
      print("Compiling for masking of order", i,"and rng", rng)
      system("make clean > /dev/null && make ORDER="+str(i)+" RNG="+str(rng)+REDIRECT)
      print("Running tests...", end ='')
      sys.stdout.flush()
      cur = popen("./main_test").read()
      print("Writing to "+PATH+" ...")
      f.write(cur)
      res += cur
      print(cur)
      print("Done.")