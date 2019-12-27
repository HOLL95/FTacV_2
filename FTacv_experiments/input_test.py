import sys, select

print("You have ten seconds to answer!")

i, o, e = select.select( [sys.stdin], [], [], 5)
print(len(i), o, e)
if (i):
  print("You said", sys.stdin.readline().strip())
else:
  print("You said nothing!")
