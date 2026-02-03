import matplotlib.pyplot as plt

def plot_molecule(coords):
    # coords: List[Dict] with keys 'index','symbol','x','y','z'
    xs = [c['x'] for c in coords]
    ys = [c['y'] for c in coords]
    zs = [c['z'] for c in coords]
    labels = [f"{c['symbol']}{c['index']}" for c in coords]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs, ys, zs)

    # annotate each point
    for x, y, z, label in zip(xs, ys, zs, labels):
        ax.text(x, y, z, label, size=8)

    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title('3D Structure from SMILES Coordinates')
    plt.tight_layout()
    plt.show()
    
