import numpy as np
import matplotlib.pyplot as plt

plt.style.use('seaborn-whitegrid')


def question_5():
    csv_data = np.genfromtxt('H3_modeling_100_scores.csv', delimiter=',',
                             skip_header=1, usecols=range(1, 3))
    # y is total score
    y = [val[0] for val in csv_data]
    # X is RMSD
    X = [val[1] for val in csv_data]

    plt.scatter(X, y, marker='o')
    plt.xlabel('RMSD')
    plt.ylabel('rosetta scores')
    plt.show()


if __name__ == "__main__":
    question_5()
