'''
    Cluster shape identification class.
    Intended to be imported and used outside of the script.
    STILL WORK IN PROGRESS - DO NOT USE AS IT IS!
    
    Author: Giorgio Alberto Lucia
    contact: giorgio.lucia@edu.unito.it
'''

import numpy as np
import cv2
from itertools import product
import matplotlib.pyplot as plt

class PatternIdentifier:
    def __init__(self, n, m, N, output_path=None):
        '''
        Initialize the PatternIdentifier class.

        Parameters
        ----------
            - n (int): Height of the patterns.
            - m (int): Width of the patterns.
            - N (int): Maximum number of ones in a pattern.
            - output_path (str, optional): Output file path to save the patterns image. Defaults to None. (WIP, do not pass this variable)
        '''
        self.n = n
        self.m = m
        self.N = N
        self.patterns = self.generatePatterns()
        self.output_path = output_path

        if self.output_path:
            self.savePatternsImage()

    def generatePatterns(self):
        '''
        Generate all possible patterns with up to N ones.

        Returns
        -------
            - list: List of generated patterns.
        '''
        patterns = []
        for num_ones in range(1, self.N + 1):
            for positions in product(range(self.n * self.m), repeat=num_ones):
                pattern = np.zeros((self.n, self.m), dtype=np.uint8)
                for pos in positions:
                    row, col = divmod(pos, self.m)
                    pattern[row, col] = 1
                patterns.append(pattern)
        return patterns

    def identifyPattern(self, input_matrix):
        '''
        Identify a pattern in an input matrix.

        Parameters
        ----------
            - input_matrix (list): Input matrix containing 0s and 1s.

        Returns
        -------
            - int: Index of the identified pattern, or -1 if not found.
        '''

        input_matrix = np.array(input_matrix, dtype=np.uint8)
        input_moments = cv2.moments(input_matrix)
        input_hu_moments = cv2.HuMoments(input_moments).flatten()

        for idx, pattern in enumerate(self.patterns):
            for rotation in range(4):
                rotated_pattern = np.rot90(pattern, rotation)
                rotated_moments = cv2.moments(rotated_pattern)
                rotated_hu_moments = cv2.HuMoments(rotated_moments).flatten()

                if np.allclose(input_hu_moments, rotated_hu_moments, rtol=1e-5, atol=1e-5):
                    return idx
        return -1
    
    def savePatternsImage(self):
        '''
        Save an image displaying all generated patterns with their indexes.
        (DOES NOT WORK)
        '''
        num_patterns = len(self.patterns)
        fig, axes = plt.subplots(1, num_patterns, figsize=(12, 3))
        
        for idx, pattern in enumerate(self.patterns):
            axes[idx].imshow(pattern, cmap='gray')
            axes[idx].set_title(f"Pattern {idx}")
            axes[idx].axis('off')
        
        plt.tight_layout()

        if self.output_path.endswith('.png'):
            plt.savefig(self.output_path, dpi=300, bbox_inches='tight')
        elif self.output_path.endswith('.pdf'):
            plt.savefig(self.output_path, format='pdf', bbox_inches='tight')
        plt.close()


if __name__ == '__main__':
    # Example usage
    n = 4
    m = 4
    N = 3
    output_path = 'patterns.png'
    
    input_matrix = [
        [1, 0, 0, 0],
        [1, 0, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 0, 0],
    ]
    input_matrix2 = [
        [1, 1, 1, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]
    input_matrix3 = [
        [0, 0, 0, 0],
        [0, 1, 1, 1],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]
    
    pattern_identifier = PatternIdentifier(n, m, N, output_path)
    pattern_id = pattern_identifier.identifyPattern(input_matrix)
    print(f"Identified Pattern ID: {pattern_id}")
    pattern_id2 = pattern_identifier.identifyPattern(input_matrix2)
    print(f"Identified Pattern ID 2: {pattern_id2}")
    pattern_id3 = pattern_identifier.identifyPattern(input_matrix3)
    print(f"Identified Pattern ID 3: {pattern_id3}")
