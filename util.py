def print_selection_of(_list):
    index = 0
    for element in _list:
        print(str(index) + ") ", element)
        index += 1


read_int = lambda msg: int(input(msg))
