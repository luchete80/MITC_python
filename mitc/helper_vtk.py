def write_vtk(filename, nodes, elements, displacements):
    """
    Write a VTK file with nodes, elements, and displacement data.

    Args:
    - filename: Name of the VTK output file.
    - nodes: List of [x, y, z] coordinates for each node.
    - elements: List of [n1, n2, n3, n4] (0-based indices) for each quad element.
    - displacements: List of [ux, uy, uz] for each node.
    """
    with open(filename, "w") as f:
        # Header
        f.write("# vtk DataFile Version 2.0\n")
        f.write("VTK Unstructured Grid Example\n")
        f.write("ASCII\n\n")

        # Points (Nodes)
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(nodes)} float\n")
        for node in nodes:
            f.write(f"{node[0]} {node[1]} {node[2]}\n")

        # Elements (Cells)
        num_elements = len(elements)
        cell_size = sum([len(e) + 1 for e in elements])  # Extra 1 for number of nodes per element
        f.write(f"\nCELLS {num_elements} {cell_size}\n")
        for element in elements:
            f.write(f"{len(element)} " + " ".join(map(str, element)) + "\n")

        # Cell types (5 = VTK_TRIANGLE for triangle elements)
        f.write(f"\nCELL_TYPES {num_elements}\n")
        for _ in elements:
            f.write("5\n")  # 5 = VTK_TRIANGLE

        # Node-based displacement data
        f.write(f"\nPOINT_DATA {len(nodes)}\n")
        f.write("VECTORS Displacement float\n")
        for disp in displacements:
            f.write(f"{disp[0]} {disp[1]} {disp[2]}\n")

        f.write("VECTORS Rotations float\n")
        for disp in displacements:
            f.write(f"{disp[3]} {disp[4]} {disp[5]}\n")

# Example Data
# nodes = [
#    [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]  # 4 nodes
# ]

# elements = [
#    [0, 1, 2, 3]  # 1 quadrilateral shell element
# ]

# displacements = [
#    [0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.1, 0.1, 0.0], [0.0, 0.1, 0.0]
# ]

# Write VTK file
# write_vtk("shell_mesh.vtk", nodes, elements, displacements)

