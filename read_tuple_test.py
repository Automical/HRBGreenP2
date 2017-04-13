#read tuple from txt test
import ast
f = open('h.txt', 'r')

l = f.readlines()[0]
print(l)
l=l.replace('(','[')
l=l.replace(')',']')
l=l.replace(' ','')
print(l)
print("")
z=ast.literal_eval(l)

print("converted")
print("")

print(z)
print("first set")
print(z[0])
print("first first element")
print(z[0][0])

print("")

print("second set")
print(z[1])
print("second first element")
print(z[1][0])

print(z.length)
print(z[0].length)

#print(z[0])
#print input(l)
