import random

def scoring_line(acc1, acc2, sc):
    print(f'{acc1}	{acc2}	0	0	0	0	0	0	0	0	0.00001	{sc}')	


def one_correlated_pair(i, j, subjects):
    for k in subjects:
        score = random.uniform(0, 200)
        scoring_line(f's{i}', f't{k}', score)
        scoring_line(f's{j}', f't{k}', score)

def one_uncorrelated_pair(i, j, subjects):
    for k in subjects:
        score = random.uniform(0, 200)
        scoring_line(f's{i}', f't{k}', score)
        score = random.uniform(0, 200)
        scoring_line(f's{j}', f't{k}', score)
        

def main():
    one_correlated_pair(0, 1, range(10))
    one_uncorrelated_pair(0, 2, range(10,15))
    one_uncorrelated_pair(1, 2, range(15,20))

if __name__ == '__main__':
    main()
        
