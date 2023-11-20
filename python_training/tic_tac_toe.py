import pprint

theBoard = {'top-L': 'O',
            'top-M': 'O',
            'top-R': 'O',
            'mid-L': 'X',
            'mid-M': 'X',
            'mid-R': 'X',
            'low-L': ' ',
            'low-M': ' ',
            'low-R': ' '}

pprint.pprint(theBoard)

def printBoard(board):
    print(board['top-L'] + '|' + board['top-M'] + '|' + board['top-R'])
    print('-----')
    print(board['mid-L'] + '|' + board['mid-M'] + '|' + board['mid-R'])
    print('-----')
    print(board['low-L'] + '|' + board['low-M'] + '|' + board['low-R'])

printBoard(theBoard)
