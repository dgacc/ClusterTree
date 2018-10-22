package operator;


import java.util.ArrayList;
import java.util.Collections;

import dislay.Paint;
import dislay.Windows;
import filesinout.ReadFiles;
import structures.Cycles;

public class CayleyCode {
	public static void main(String[] args) {
		ReadFiles.clusterReadFiles("5eil76.clt");
		Windows windows = new Windows();
		windows.runWindow(" cac cay");
		Paint p = new Paint();
		CayleyCode cc = new CayleyCode();
		// 17, 7, 7, 3, 12, 18, 1, 8, 20, 1, 6, 13,
		// 4, 17, 4, 13, 5, 10
		// int[] dandelionCode = { 17, 7, 7, 3, 12, 18, 1, 8, 19, 1, 6, 13, 4,
		// 17, 4, 13, 5, 10 };
		// int[] dandelionCode = { 16, 6, 6, 2, 11, 17, 0, 7, 19, 0, 5, 12, 3,
		// 16, 3, 12, 4, 9 };
		// 17, 5, 7, 3, 13, 1, 8, 1, 20, 4, 6, 18, 4, 17,7, 13, 12, 10
		int[] blobCode = { 16, 4, 6, 2, 12, 0, 7, 0, 19, 3, 5, 17, 3, 16, 6, 12, 11, 9 };
		// int[] blobCode = { 4, 6, 12, 0, 7, 0, 3, 5, 3, 6, 12, 11 }; // int[]
		// dandelionCode
		// = {5, 3,
		// 3, 6, 9,
		// 5, 6, 3,
		// 3, 7};
		double[][] caytest = cc.decodingBlob(blobCode);

		p.setPaint(caytest, ReadFiles.vertices, ReadFiles.clusters, 20, 0, 0, ReadFiles.root);
		windows.addPaint(p);

	}

	public void endcodingDandelion(double[][] Tree) {

	}

	/**
	 * In order to construct the tree T corresponding to DanelionCode D
	 * 
	 * @param DandelionCode
	 *            a Dandelion String D
	 * @return The tree corresponding to D
	 */
	public double[][] decodingDandelion(int[] DandelionCode) {

		// define the mapping function
		int DandelionLength = DandelionCode.length;
		int[] tempTable = new int[DandelionLength];
		boolean[] mask = new boolean[DandelionLength + 2];
		double[][] tree = new double[DandelionLength + 2][DandelionLength + 2];

		for (int i = 0; i < DandelionLength; i++) {
			tempTable[i] = i + 1;
			mask[i] = false;
		}
		mask[DandelionLength + 1] = false;
		mask[DandelionLength] = false;
		// step 2 + 3: Let the cycles associated with the function D.
		// determine cycles;

		ArrayList<Cycles> cycle1 = new ArrayList<Cycles>();
		ArrayList<Integer> temp = new ArrayList<>();
		for (int i = 0; i < DandelionLength; i++) {
			boolean flag = false;
			int maxElement = tempTable[i];
			// Cycles cc = new Cycles();
			Cycles cycle = new Cycles();
			cycle.getCycle().add(tempTable[i]);

			int k = i;
			if (DandelionCode[k] == 0 || DandelionCode[k] == DandelionLength + 1) {
				mask[tempTable[i]] = true;
				continue;
			}
			while (!mask[DandelionCode[k]] && !temp.contains(DandelionCode[k])) {
				if (DandelionCode[k] == 0 || DandelionCode[k] == DandelionLength + 1) {
					break;
				}
				if (DandelionCode[k] > maxElement) {
					maxElement = DandelionCode[k];
				}

				if (DandelionCode[k] == cycle.getCycle().get(0)) {
					flag = true;
					break;

				} else {
					cycle.getCycle().add(DandelionCode[k]);
					mask[DandelionCode[k]] = true;
					k = DandelionCode[k] - 1;

				}
			}
			mask[tempTable[i]] = true;
			// is this satisfies condition, if it is not satisfy conditions then
			// return to the the first state
			if (!flag) {
				for (int j = 0; j < cycle.getCycle().size(); j++) {
					mask[cycle.getCycle().get(j)] = false;
				}
			} else {
				// sort element from the cycle
				int index = cycle.getCycle().indexOf(maxElement);
				int length = cycle.getCycle().size();

				for (int j = length - 1; j > index; j--) {
					if (cycle.getCycle().get(j) < maxElement) {
						cycle.getCycle().addFirst(cycle.getCycle().get(j));
						cycle.getCycle().remove(length);
					}
				}
				cycle.setMaxElement(maxElement);
				cycle1.add(cycle);
				temp.addAll(cycle.getCycle());

			}
		}
		Collections.sort(cycle1, ChromosomeCmp.compareByMaxElement);
		ArrayList<Integer> path = new ArrayList<>();
		int numberOfCycles = cycle1.size();
		if (numberOfCycles > 0) {
			for (int i = 0; i < numberOfCycles; i++) {
				path.addAll(cycle1.get(i).getCycle());
			}
			int piLength = path.size();
			// step 4: Construct the Tree corresponding to D
			// 4.1: create a path from vertex 1 to vertex n by following the
			// list in
			// pi from left to right
			tree[0][path.get(0)] = tree[path.get(0)][0] = 1.0f;
			tree[DandelionLength + 1][path
					.get(piLength - 1)] = tree[path.get(piLength - 1)][DandelionLength + 1] = 1.0f;
			for (int i = 0; i < piLength - 1; i++) {
				tree[path.get(i)][path.get(i + 1)] = tree[path.get(i + 1)][path.get(i)] = 1.0f;
			}
		} else {
			tree[0][DandelionLength + 1] = tree[DandelionLength + 1][0] = 1.0f;
		}
		// 4.2 the others vertex are not in pi
		for (int i = 0; i < DandelionLength; i++) {
			if (!path.contains(i + 1)) {
				tree[i + 1][DandelionCode[i]] = tree[DandelionCode[i]][i + 1] = 1.0f;
			}
		}
		return tree;

	}

	public void endcodingBlob(double[][] Tree, int numVertices) {
		// step 1: let all of the vertices be uncolored set v = n
		int v = numVertices - 1;
		int[] state = new int[numVertices];// values meant 0: uncolored, 1:
											// black, 2: white
		// step 2: Assign the color back to vertex v
		state[v] = 1;
		// step 3: Assign the color white to every uncolored ancestor of v in
		// the tree T
		for (int i = 0; i < numVertices; i++) {
			// Assign the color white;
		}
		// step 4: if

	}

	/**
	 * In order to Construct the Tree corresponding to the blobcode
	 * 
	 * @param Blobcode
	 *            is a Blob String B
	 * @return The tree corresponding to BlobCode B
	 */
	public double[][] decodingBlob(int[] Blobcode) {

		int BlobcodeLength = Blobcode.length;
		int numVertex = BlobcodeLength + 2;
		double[][] tree = new double[numVertex][numVertex];
		int[] color = new int[numVertex];// values meant 0: uncolored,
		// step 2 + 3: Let the cycles associated with the function D.
		// determine cycles;
		int k = 0;
		for (int i = 0; i < BlobcodeLength; i++) {
			// tim dinh trang:
			k = i;
			boolean[] mask = new boolean[numVertex];
			mask[i + 1] = true;
			boolean flag = false;
			if (Blobcode[i] == 0) {
				color[i + 1] = 1;
				continue;
			}
			int vertex = i + 1;
			while ((vertex) >= Blobcode[k]) {
				mask[Blobcode[k]] = true;
				if (Blobcode[k] == 0) {
					flag = true;
					color[i + 1] = 1;
					break;
				}
				k = Blobcode[k] - 1;
				if (mask[Blobcode[k]]) {
					color[i + 1] = 1;
					flag = true;
					break;
				}
			}
			if (!flag) {
				color[i + 1] = 2;
			}
		}
		color[0] = 1;
		color[BlobcodeLength + 1] = 1;

		// Step 5: Suppose the black vertices are x1 < x2 < . . . < xt,
		// where t ∈ [2, n] is the total number of black vertices; observe
		// that x1 = 1 and xt = n. To construct the tree T ∈ Tn
		// corresponding to B, take a set of n isolated vertices (labelled
		// with the integers from 1 to n), create the edge (i, bi) for each
		// white vertex i ∈ [2, n − 1], create the edge (xi, bxi−1) for
		// each i ∈ [3, t], and finally create the edge (x2, 1)
		ArrayList<Integer> blackVertices = new ArrayList<>();
		for (int i = 0; i < numVertex; i++) {
			if (color[i] == 1) {
				blackVertices.add(i);
			}
		}
		// blackVertices.add(0);
		Collections.sort(blackVertices);
		// create the edge (i, bi) for each white vertex i ∈ [2, n − 1]
		for (int i = 1; i < BlobcodeLength + 1; i++) {
			if (color[i] == 2) {
				tree[i][Blobcode[i - 1]] = 1.0f;
				tree[Blobcode[i - 1]][i] = 1.0f;

			}
		}

		// create the edge (xi, bxi−1) for each i ∈ [3, t]
		int t = blackVertices.size();
		int xi = 0;
		int xi_1 = 0;
		for (int i = 2; i < t; i++) {
			xi = blackVertices.get(i);
			xi_1 = blackVertices.get(i - 1);
			tree[xi][Blobcode[xi_1 - 1]] = 1.0f;
			tree[Blobcode[xi_1 - 1]][xi] = 1.0f;
		}
		// and finally create the edge (x2, 1)
		tree[blackVertices.get(1)][0] = 1.0f;
		tree[0][blackVertices.get(1)] = 1.0f;
		return tree;
	}

	public boolean AreAllColored(int[] color, int colorLength) {
		boolean flag = true;
		for (int i = 0; i < colorLength; i++) {
			if (color[i] == 0) {
				flag = false;
				break;
			}
		}
		return flag;
	}
}
