class RunModifier:

    def __init__(self, name, modifier_function, count=5):
        self.name = name
        self.index = count
        self.mod = modifier_function

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == 0:
            raise StopIteration
        self.index = self.index-1
        return self.index

    def modifier(self):
        return self.modifier_function(self.index)