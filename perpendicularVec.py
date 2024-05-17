import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def perpendicular_vectors(v, n):
    # Ensure the input vector is normalized
    v = v / np.linalg.norm(v)
    
    # Generate an arbitrary vector not parallel to v
    not_parallel = np.array([1.0, 0.0, 0.0])
    if np.dot(v, not_parallel) == 1.0:
        not_parallel = np.array([0.0, 1.0, 0.0])
        
    #the first perpendicular vector by taking the cross product
    v1 = np.cross(v, not_parallel)
    v1 = v1 / np.linalg.norm(v1)
    
    #the second perpendicular vector by taking the cross product again
    v2 = np.cross(v, v1)
    v2 = v2 / np.linalg.norm(v2)
    
    # Initialize a list to store the perpendicular vectors
    perpendicular_vectors = [v1, v2]
    
    # For additional perpendicular vectors, rotate v1 and v2 by 360/n degrees around v
    for i in range(2, n * 2, 2):
        angle = i * (360 / (n * 2))
        rotation_matrix = rotation_matrix_from_vectors(v, angle)
        rotated_v1 = np.dot(rotation_matrix, v1)
        rotated_v2 = np.dot(rotation_matrix, v2)
        perpendicular_vectors.append(rotated_v1)
        perpendicular_vectors.append(rotated_v2)
        
    return perpendicular_vectors

def rotation_matrix_from_vectors(axis, theta):
    """
    Compute the rotation matrix to rotate around an arbitrary axis by a given angle.
    """
    axis = axis / np.linalg.norm(axis)
    a = np.cos(np.deg2rad(theta) / 2.0)
    b, c, d = -axis * np.sin(np.deg2rad(theta) / 2.0)
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])


if __name__=="__main__":
    # Example usage:
    v = np.array([1, 2, 0])  # Example input vector
    n = 15  # Number of unit vectors perpendicular to v to generate
    perpendicular_vectors_list = perpendicular_vectors(v, n)

    # Plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Original vector
    ax.quiver(0, 0, 0, v[0], v[1], v[2], color='blue', label='Original Vector')

    # Perpendicular vectors
    for i, vec in enumerate(perpendicular_vectors_list):
        ax.quiver(0, 0, 0, vec[0], vec[1], vec[2], color='red', label=f'Perpendicular Vector {i+1}')

    # Set plot limits
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])

    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Perpendicular Vectors')

    # Add legend
    ax.legend()

    plt.show()
