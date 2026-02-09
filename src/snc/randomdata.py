import random

def scoring_line(acc1, acc2, sc):
    print(f'{acc1}	{acc2}	0	0	0	0	0	0	0	0	0.00001	{sc}')	


def one_correlated_pair(i, j, subjects):
    for k in subjects:
        score = random.uniform(0, 200)
        scoring_line(f's{i}', f't{k}', score)
        scoring_line(f's{j}', f't{k}', score)

def one_uncorrelated(i, subjects):
    for k in subjects:
        score = random.uniform(0, 200)
        scoring_line(f's{i}', f't{k}', score)
        

def main():
    one_correlated_pair(0, 1, range(10))
    one_uncorrelated(2, range(10))

if __name__ == '__main__':
    main()
        
