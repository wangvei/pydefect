#!/usr/bin/python

class Hoge:
    def __init__(self, val):
        self.val = val

    def change(self):
        self.val += 1

    def __repr__(self):
        return str(self.val)

l = []
for i in range(10):
    l.append(Hoge(i))
print(l)
for hoge in l:
    hoge.change()
print(l)
