# change 0.0 to 0
from matplotlib.ticker import FuncFormatter


# need pos as the FuncFormatter arg requires 2 args.
# https://stackoverflow.com/questions/29139552/pyplot-remove-the-digits-of-zero-start-from-0-not-0-00
def my_formatter(number, pos):
    return str(int(number)) if number.is_integer() else str(number)


formatter = FuncFormatter(my_formatter)
