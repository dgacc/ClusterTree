package operator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;
import filesinout.ReadFiles;
import pruferpresentation.MFEAPruferNumber;
import random.MyRandom;
import structures.Cluster;
import structures.Edge;
import structures.Individual;

public class InitializeChromosome {
	public static ArrayList<Cluster> maxClusters = new ArrayList<Cluster>();
	// public static boolean isInstance1[]; // is number of max vertices in
	// cluster equals number of vertex in cluster instance1
	// public static boolean isMaxClusterInstance1; //Return true if cluster in
	// instance 1 have the number of vertex max.
	// private Random a = new Random();
	static int minNumberOfCluster = 0; // number of cluster of less between two
										// instances
	static int maxNumberOfCluster = 0; // number of cluster of greater between
										// two instances
	public static int maxNumberVertices = 0; // number of vertices of
												// maxInstances
	public static int[] maxClusterVertices; // array of number of vertex in
											// cluster have greater than the
											// others.
	public static int[] minClusterVertices; // array of number of vertex in
											// cluster have less than the
											// others.

	/**
	 * copy the weightMatrix for each Cluster
	 * 
	 * @param weightMatrix
	 * @param cluster
	 * @return cluster weightMatrix
	 */
	public double[][] buildClusterWeightMatrix(double[][] weightMatrix, ArrayList<Integer> cluster) {
		int numberClusterVertices = cluster.size();
		double[][] clusterWeightMatrix = new double[numberClusterVertices][numberClusterVertices];
		for (int i = 0; i < numberClusterVertices; i++) {
			for (int j = 0; j < numberClusterVertices; j++) {
				clusterWeightMatrix[i][j] = weightMatrix[cluster.get(i)][cluster.get(j)];
			}
		}
		return clusterWeightMatrix;
	}

	/**
	 * generate randomly tree based on PRimRST first: generate the tree for each
	 * cluster then convert to the TRee Second: generate the presentation tree:
	 * meaning each cluster is presented by a vertex then covert to the Tree
	 * 
	 * @param weightMatrix
	 * @param clusters
	 * @param num_vertices
	 * @return Tree after do that stuff
	 */
	public double[][] clusterPrimRST(double[][] weightMatrix, ArrayList<Cluster> clusters, int num_vertices) {
		// int randomVertice = -1, indexRandomVertice = -1;
		double[][] tree = new double[num_vertices][num_vertices];
		double[][] clusterWeightMatrix;
		double[][] clusterSpanningTree;
		// int[] indexPresentationVertex = new int[clusters.size()];
		// build the spanning for each cluster
		for (int i = 0; i < ReadFiles.numberOfCluster; i++) {

			int numberClusterVertices = clusters.get(i).getCluster().size();
			clusterWeightMatrix = buildClusterWeightMatrix(weightMatrix, clusters.get(i).getCluster());
			clusterSpanningTree = primRST(clusterWeightMatrix, numberClusterVertices);

			for (int j = 0; j < numberClusterVertices; j++) {
				for (int k = 0; k < numberClusterVertices; k++) {
					if (clusterSpanningTree[j][k] > 0) {
						tree[clusters.get(i).getCluster().get(j)][clusters.get(i).getCluster()
								.get(k)] = clusterSpanningTree[j][k];
					}
				}
			}
		}

		// Each cluster is represented by one vertex
		Cluster clusterRepresentation = new Cluster();
		for (int i = 0; i < ReadFiles.numberOfCluster; i++) {
			int vertex = MyRandom.r.nextInt(clusters.get(i).getCluster().size());
			clusterRepresentation.getCluster().add(clusters.get(i).getCluster().get(vertex));
		}
		clusterWeightMatrix = buildClusterWeightMatrix(weightMatrix, clusterRepresentation.getCluster());
		clusterSpanningTree = primRST(clusterWeightMatrix, ReadFiles.numberOfCluster);
		// Transform to spanning tree of G Graph
		for (int i = 0; i < ReadFiles.numberOfCluster; i++) {
			for (int j = 0; j < ReadFiles.numberOfCluster; j++) {
				if (clusterSpanningTree[i][j] > 0) {
					tree[clusterRepresentation.getCluster().get(i)][clusterRepresentation.getCluster()
							.get(j)] = clusterSpanningTree[i][j];
				}
			}
		}
		return tree;
	}

	/**
	 * this is algorithms to generate the tree randomly first: apply "findEdge"
	 * method to find out the list of addjence
	 * 
	 * @param weightMatrix
	 * @param num_vertices
	 * @return
	 */
	public double[][] primRST(double[][] weightMatrix, int num_vertices) {
		int randomVertice = -1, indexRandomEdge = -1;
		ArrayList<Integer> contain = new ArrayList<Integer>();
		ArrayList<Edge> edgeListAdjacence = new ArrayList<Edge>();
		double[][] tempTable = new double[num_vertices][num_vertices];

		randomVertice = MyRandom.r.nextInt(num_vertices);
		contain.add(randomVertice);
		edgeListAdjacence.addAll(findEdge(num_vertices, randomVertice, contain, weightMatrix));

		while (contain.size() < num_vertices) {

			if (contain.size() == 0) {
				return null;
			}

			indexRandomEdge = MyRandom.r.nextInt(edgeListAdjacence.size());
			Edge tempEdge = new Edge(edgeListAdjacence.get(indexRandomEdge).startVertice,
					edgeListAdjacence.get(indexRandomEdge).endVertice);
			edgeListAdjacence.remove(indexRandomEdge);

			if (!contain.contains(tempEdge.endVertice)) {
				//
				tempTable[tempEdge.startVertice][tempEdge.endVertice] = 1.0f;
				tempTable[tempEdge.endVertice][tempEdge.startVertice] = 1.0f;
				contain.add(tempEdge.endVertice);
				edgeListAdjacence.addAll(findEdge(num_vertices, tempEdge.endVertice, contain, weightMatrix));
			}
		}
		return tempTable;
	}

	/**
	 * find that addjence with that edge
	 * 
	 * @param num_vertices
	 * @param vertice
	 * @param contain
	 * @param weightMatrix
	 * @return the list of the addjence edge
	 */
	public ArrayList<Edge> findEdge(int num_vertices, int vertice, ArrayList<Integer> contain,
			double[][] weightMatrix) {
		ArrayList<Edge> listEdge = new ArrayList<Edge>();
		for (int i = 0; i < num_vertices; i++) {
			if (i != vertice) {
				if (!contain.contains(i)) {
					if (weightMatrix[vertice][i] > 0) {
						Edge tempEdge = new Edge(vertice, i);
						listEdge.add(tempEdge);

					}
				}
			}
		}
		return listEdge;
	}

	// get a Edge from the list of edges
	public Edge getEdge(int index, ArrayList<Edge> edgeList) {
		return edgeList.get(index);
	}

	/**
	 * 
	 * @param start_vertice
	 * @param end_vertice
	 * @param number_vertices
	 * @param weight_matrix
	 * @param visited
	 * @param pre
	 */
	public void findCycle(int start_vertice, int end_vertice, int number_vertices, double[][] weight_matrix,
			boolean[] visited, int[] pre) {
		visited[start_vertice] = true;
		for (int i = 0; i < number_vertices; i++) {
			if ((weight_matrix[start_vertice][i] > 0) && (i == end_vertice)) {
				pre[i] = start_vertice;
				return;
			}
			if ((weight_matrix[start_vertice][i] > 0) && (!visited[i])) {
				pre[i] = start_vertice;
				findCycle(i, end_vertice, number_vertices, weight_matrix, visited, pre);
			}
		}
	}

	/**
	 * find the tree based on BFS, make a tree
	 * 
	 * @param weightMatrix
	 * @param num_vertex
	 * @param startVertex
	 * @return BFS Tree
	 */
	public double[][] BFSTree(double[][] weightMatrix, int num_vertex, int startVertex) {

		double[][] BFStree = new double[num_vertex][num_vertex];
		boolean[] mark = new boolean[num_vertex];
		Queue<Integer> queue = new LinkedList<>();
		// initialize label for each vertex,
		for (int i = 0; i < num_vertex; i++) {
			mark[i] = true;
		}
		mark[startVertex] = false;
		queue.add(startVertex);

		while (!queue.isEmpty()) {
			int considerVertex = queue.poll();
			for (int i = 0; i < num_vertex; i++) {
				if (weightMatrix[considerVertex][i] > 0 && mark[i]) {
					queue.add(i);
					mark[i] = false;
					BFStree[considerVertex][i] = 1.0f;
					BFStree[i][considerVertex] = 1.0f;

				}
			}

		}
		return BFStree;
	}

	/**
	 * can sua dau vao cho bai toan
	 * 
	 * @param numCluster
	 * @param weightMatrix
	 * @param vertex_In_Cluster
	 * @param idx
	 * @param edMat
	 * @return
	 */

	public double[][] localSearch(int numCluster, double[][] weightMatrix, int[][] vertex_In_Cluster, int[] idx,
			double[][] edMat) {
		double temp[][] = new double[weightMatrix.length][weightMatrix.length];
		for (int i = 0; i < weightMatrix.length; i++) {
			for (int j = 0; j < weightMatrix.length; j++) {
				temp[i][j] = edMat[i][j];
			}
		}
		for (int i = 0; i < numCluster; i++) {
			for (int j = i + 1; j < numCluster; j++) {
				if (edMat[idx[i]][idx[j]] > 0) {
					int[] arr = findMinDistance(weightMatrix, vertex_In_Cluster[i], vertex_In_Cluster[j]);
					temp[idx[i]][idx[j]] = 0;
					idx[i] = arr[0];
					idx[j] = arr[1];
					temp[idx[i]][idx[j]] = 1;
				}
			}
		}
		return temp;
	}

	public int[] findMinDistance(double[][] weightMatrix, int[] a1, int a2[]) {
		int result[] = new int[2];
		double min = weightMatrix[a1[0]][a2[0]];
		for (int i = 0; i < a1.length; i++) {
			for (int j = 0; j < a2.length; j++) {
				if (weightMatrix[a1[i]][a2[j]] < min) {
					min = weightMatrix[a1[i]][a2[j]];
					result[0] = a1[i];
					result[1] = a2[j];
				}
			}
		}
		return result;
	}

	/**
	 * After call this method all of the Information in the common gene going to
	 * be updated Calculate the number of vertices in cluster min for each
	 * cluster( Calculate the number of vertices in cluster max for each cluster
	 * And calculate the maxNumberVertices, then calculate the number of cluster
	 * for each instance Update the value of each the instance variables on top
	 * 
	 * @param clusters1
	 *            : The list of cluster in instance 1
	 * @param clusters2
	 *            : The list of cluster in instance 2
	 */
	public void clutersVerticesInf(ArrayList<Cluster> clusters1, ArrayList<Cluster> clusters2) {
		maxNumberVertices = 0;
		if (clusters1.size() <= clusters2.size()) {
			minNumberOfCluster = clusters1.size();
			maxNumberOfCluster = clusters2.size();
			minClusterVertices = new int[maxNumberOfCluster]; // the number of
																// vertices less
																// than in two
																// instances
			maxClusterVertices = new int[maxNumberOfCluster]; // the number of
																// vertices
																// greater than
																// in two
																// instances

			for (int i = 0; i < minNumberOfCluster; i++) {
				if (clusters1.get(i).getCluster().size() < clusters2.get(i).getCluster().size()) {
					minClusterVertices[i] = clusters1.get(i).getCluster().size();
					maxClusterVertices[i] = clusters2.get(i).getCluster().size();
					maxNumberVertices += maxClusterVertices[i];
				} else {
					minClusterVertices[i] = clusters2.get(i).getCluster().size();
					maxClusterVertices[i] = clusters1.get(i).getCluster().size();
					maxNumberVertices += maxClusterVertices[i];
				}
			}

			for (int j = minNumberOfCluster; j < clusters2.size(); j++) {
				minClusterVertices[j] = maxClusterVertices[j] = clusters2.get(j).getCluster().size();
				maxNumberVertices += maxClusterVertices[j];
			}
		} else {
			minNumberOfCluster = clusters2.size();
			maxNumberOfCluster = clusters1.size();

			for (int i = 0; i < minNumberOfCluster; i++) {
				if (clusters1.get(i).getCluster().size() < clusters2.get(i).getCluster().size()) {
					minClusterVertices[i] = clusters1.get(i).getCluster().size();
					maxClusterVertices[i] = clusters2.get(i).getCluster().size();
					maxNumberVertices += maxClusterVertices[i];
				} else {
					minClusterVertices[i] = clusters2.get(i).getCluster().size();
					maxClusterVertices[i] = clusters1.get(i).getCluster().size();
					maxNumberVertices += maxClusterVertices[i];
				}
			}

			for (int j = minNumberOfCluster; j < clusters1.size(); j++) {
				minClusterVertices[j] = maxClusterVertices[j] = clusters1.get(j).getCluster().size();
				maxNumberVertices += maxClusterVertices[j];
			}
		}

	}

	/**
	 * 
	 * @param maxClusters
	 *            : the cluster is mapped to use for two instances( List of
	 *            clusters in two instances)
	 * @param r
	 * @return
	 */

	public double[][] initializeGroupGene(ArrayList<Cluster> maxClusters, Random r) {

		double[][] treeMax = new double[maxNumberVertices][maxNumberVertices];
		double[][] weightMatrix = new double[maxNumberVertices][maxNumberVertices];
		for (int i = 0; i < maxNumberVertices; i++) {
			for (int j = 0; j < maxNumberVertices; j++) {
				weightMatrix[i][j] = 1.0f;
			}
		}

		// initialize gene for each cluster
		for (int i = 0; i < maxNumberOfCluster; i++) {

			double[][] bigClusterTree = new double[maxClusterVertices[i]][maxClusterVertices[i]];
			double[][] smallClusterTree = new double[minClusterVertices[i]][minClusterVertices[i]];

			for (int j = 0; j < minClusterVertices[i]; j++) {
				for (int k = 0; k < minClusterVertices[i]; k++) {
					smallClusterTree[j][k] = 1.0f;
				}
			}
			// use primRST to generate tree for smallClusterTree
			smallClusterTree = primRST(smallClusterTree, minClusterVertices[i]);
			// copy smallClusterTree to bigClusterTree
			for (int j = 0; j < minClusterVertices[i]; j++) {
				for (int k = 0; k < minClusterVertices[i]; k++) {
					bigClusterTree[j][k] = smallClusterTree[j][k];
				}
			}

			// add edges to bigClusterTree without making cycle randomly.
			int tempNum = minClusterVertices[i];
			while (tempNum < maxClusterVertices[i]) {
				int randomIndex = r.nextInt(tempNum);
				bigClusterTree[randomIndex][tempNum] = 1.0f;
				bigClusterTree[tempNum][randomIndex] = 1.0f;
				tempNum++;
			}

			// transform to the max tree( Tree of Two instances)
			for (int j = 0; j < maxClusterVertices[i]; j++) {
				for (int k = 0; k < maxClusterVertices[i]; k++) {
					if (bigClusterTree[j][k] > 0) {
						treeMax[maxClusters.get(i).getCluster().get(j)][maxClusters.get(i).getCluster()
								.get(k)] = bigClusterTree[j][k];
					}
				}
			}

		}

		// ---------------------------End initialize gene for each
		// cluster---------------------------//

		// -----------------------Start to initialize for presentation
		// cluster----------------------//

		// build up the Tree for presentation Cluster in instance which has the
		// number of clusters less than.
		double[][] clusterWeightMatrix;
		double[][] smallClusterSpanningTree;
		// choose randomly the presentation vertex is in instance have less than
		// number of clusters
		Cluster clusterRepresentation = new Cluster();
		for (int i = 0; i < minNumberOfCluster; i++) {
			int vertex = r.nextInt(minClusterVertices[i]);
			clusterRepresentation.getCluster().add(maxClusters.get(i).getCluster().get(vertex));

		}

		clusterWeightMatrix = buildClusterWeightMatrix(weightMatrix, clusterRepresentation.getCluster());
		smallClusterSpanningTree = primRST(clusterWeightMatrix, minNumberOfCluster);

		// choose randomly the presentation vertex is in instance have greater
		// number of clusters
		for (int i = minNumberOfCluster; i < maxNumberOfCluster; i++) {
			int vertex = r.nextInt(minClusterVertices[i]);
			clusterRepresentation.getCluster().add(maxClusters.get(i).getCluster().get(vertex));

		}

		// copy
		double[][] bigClusterSpanningTree = new double[maxNumberOfCluster][maxNumberOfCluster];
		for (int j = 0; j < minNumberOfCluster; j++) {
			for (int k = 0; k < minNumberOfCluster; k++) {
				bigClusterSpanningTree[j][k] = smallClusterSpanningTree[j][k];
			}
		}

		// build up the Tree for presentation Cluster in instance which has the
		// number of clusters greater.
		int tempNum = minNumberOfCluster;
		while (tempNum < maxNumberOfCluster) {
			int randomIndex = r.nextInt(tempNum);
			bigClusterSpanningTree[randomIndex][tempNum] = 1.0f;
			bigClusterSpanningTree[tempNum][randomIndex] = 1.0f;
			tempNum++;
		}

		// Transform presentation cluster to TreeMax
		for (int i = 0; i < maxNumberOfCluster; i++) {
			for (int j = 0; j < maxNumberOfCluster; j++) {
				if (bigClusterSpanningTree[i][j] > 0) {
					treeMax[clusterRepresentation.getCluster().get(i)][clusterRepresentation.getCluster()
							.get(j)] = bigClusterSpanningTree[i][j];
				}
			}
		}
		return treeMax;
	}

	/**
	 * we build cluster max following by: the first cluster = (0, 1, 2 .....
	 * (number of vertices in 1st cluster max)) the second cluster = ((number of
	 * vertices in 1st cluster max), (number of vertices in 1st cluster max +
	 * 1).....) ..etc... the last cluster = (..........(number of vertices in
	 * last cluster max - 1)+(number of vertices in last cluster max ))
	 */
	public void buildCluster() {
		maxClusters.clear();
		int count = 0;

		for (int i = 0; i < maxNumberOfCluster; i++) {
			Cluster cluster = new Cluster();
			for (int j = 0; j < maxClusterVertices[i]; j++) {
				cluster.addElement(count, j);
				count++;
			}
			maxClusters.add(cluster);
		}
	}

	/**
	 * Tách từ cây biểu diễn chung ra các l�?i giải để dánh giá.
	 * 
	 * @param maxTree:
	 *            Tree of two instances
	 * @param clusters1
	 *            : List of clusters in each cluster of instance 1
	 * @param clusters2
	 *            : List of clusters in each cluster of instance 2
	 * @param num_vertices1
	 *            : Number of vertices in instance 1
	 * @param num_vertices2
	 *            : Number of vertices in instance 2
	 * @return
	 */
	public ArrayList<double[][]> decodingTwoTree(double[][] maxTree, ArrayList<Cluster> clusters1,
			ArrayList<Cluster> clusters2, int num_vertices1, int num_vertices2, Random r) {
		ArrayList<double[][]> list = new ArrayList<double[][]>();
		double[][] weightTree1 = new double[num_vertices1][num_vertices1];
		double[][] weightTree2 = new double[num_vertices2][num_vertices2];
		double[][] TreeInst1 = new double[num_vertices1][num_vertices1];
		double[][] TreeInst2 = new double[num_vertices2][num_vertices2];
		int numberOfcluster1 = clusters1.size();
		int numberOfcluster2 = clusters2.size();

		for (int i = 0; i < numberOfcluster1; i++) {
			double[][] TreeOfEachCluster = buildClusterWeightMatrix(maxTree, maxClusters.get(i).getCluster());

			// copy the tree in Each cluster in for Tree in Instance1.
			for (int j = 0; j < clusters1.get(i).getCluster().size(); j++) {
				for (int k = 0; k < clusters1.get(i).getCluster().size(); k++) {
					if (TreeOfEachCluster[j][k] > 0) {
						TreeInst1[clusters1.get(i).getCluster().get(j)][clusters1.get(i).getCluster()
								.get(k)] = TreeInst1[clusters1.get(i).getCluster().get(k)][clusters1.get(i).getCluster()
										.get(j)] = 1.0f;
					}
				}
			}
		}

		weightTree1 = findSpanningTreeBetweenClusters(maxTree, num_vertices1, maxClusters, clusters1);
		for (int j = 0; j < num_vertices1; j++) {
			for (int k = 0; k < num_vertices1; k++) {
				if (weightTree1[j][k] > 0) {
					TreeInst1[j][k] = weightTree1[j][k];
				}
			}
		}
		list.add(TreeInst1);

		for (int i = 0; i < numberOfcluster2; i++) {
			double[][] TreeOfEachCluster = buildClusterWeightMatrix(maxTree, maxClusters.get(i).getCluster());
			// copy the tree in Each cluster for Tree in Instance2.
			for (int j = 0; j < clusters2.get(i).getCluster().size(); j++) {
				for (int k = 0; k < clusters2.get(i).getCluster().size(); k++) {
					if (TreeOfEachCluster[j][k] > 0) {
						TreeInst2[clusters2.get(i).getCluster().get(j)][clusters2.get(i).getCluster()
								.get(k)] = TreeInst2[clusters2.get(i).getCluster().get(k)][clusters2.get(i).getCluster()
										.get(j)] = 1.0f;
					}
				}
			}
		}

		weightTree2 = findSpanningTreeBetweenClusters(maxTree, num_vertices2, maxClusters, clusters2);
		for (int j = 0; j < num_vertices2; j++) {
			for (int k = 0; k < num_vertices2; k++) {
				if (weightTree2[j][k] > 0) {
					TreeInst2[j][k] = weightTree2[j][k];
				}
			}
		}
		list.add(TreeInst2);

		return list;
	}

	/**
	 * find the edges connect all of clusters altogether.
	 * 
	 * @param parent
	 *            : tree input ( tree of two instances is presented by matrix
	 *            weight)
	 * @param num_vertex
	 *            : the number of vertices in tree input.
	 * @param maxClusters
	 *            : List of vertices in two instances
	 * @return Edges between clusters is separated in matrix weight
	 */
	public double[][] findSpanningTreeBetweenClusters(double[][] parent, int num_vertex, ArrayList<Cluster> maxClusters,
			ArrayList<Cluster> clusters1) {
		int numberOfCluster = clusters1.size();
		double[][] tree = new double[num_vertex][num_vertex];
		int numberClusterVertex1 = 0;
		int numberClusterVertex2 = 0;
		for (int i = 0; i < numberOfCluster; i++) {
			for (int j = 0; j < numberOfCluster; j++) {
				tree[i][j] = 0;
				tree[j][i] = 0;
			}
		}

		for (int i = 0; i < numberOfCluster; i++) {
			numberClusterVertex1 = clusters1.get(i).getCluster().size();
			for (int j = 0; j < numberClusterVertex1; j++) {
				for (int k = i + 1; k < numberOfCluster; k++) {
					numberClusterVertex2 = clusters1.get(k).getCluster().size();

					for (int t = 0; t < numberClusterVertex2; t++) {
						if (parent[maxClusters.get(i).getCluster().get(j)][maxClusters.get(k).getCluster()
								.get(t)] > 0) {

							tree[clusters1.get(i).getCluster().get(j)][clusters1.get(k).getCluster().get(t)] = 1.0f;
							tree[clusters1.get(k).getCluster().get(t)][clusters1.get(i).getCluster().get(j)] = 1.0f;
						}
					}
				}
			}
		}
		return tree;
	}

	/**
	 * 
	 * @param clusters
	 * @param popLength
	 * @param rnd
	 * @return
	 */
	public ArrayList<Individual> initiaizePopulation(ArrayList<Cluster> clusters, int popLength, Random rnd) {
		ArrayList<Individual> population = new ArrayList<Individual>();
		MFEAPruferNumber mfeaPruferNumber = new MFEAPruferNumber();
		for (int m = 0; m < popLength; m++) {
			ArrayList<Integer> temp = new ArrayList<Integer>();
			ArrayList<Integer> temp1 = new ArrayList<Integer>();
			ArrayList<Integer> temp3 = new ArrayList<Integer>();
			temp3.removeAll(temp3);
			Individual individual = new Individual();

			int numberOfClusters = clusters.size();
			int[] presentVertex = new int[numberOfClusters];

			for (int i = 0; i < numberOfClusters; i++) {

				temp.add(i, i);
				int numberClusterVertex = clusters.get(i).getCluster().size();
				presentVertex[i] = rnd.nextInt(numberClusterVertex);
				for (int j = 0; j < numberClusterVertex; j++) {
					temp1.add(j, j);
				}
				Collections.shuffle(temp1, rnd);
				temp1 = mfeaPruferNumber.deletePruferNumber(temp1, rnd);
				temp3.addAll(temp1);
				temp1.removeAll(temp1);
			}

			Collections.shuffle(temp, rnd);
			temp = mfeaPruferNumber.deletePruferNumber(temp, rnd);
			temp3.addAll(temp);
			temp.removeAll(temp);
			int num = temp3.size() + numberOfClusters;
			int k = 0;
			for (int i = 0; i < num; i++) {
				if (i < temp3.size()) {
					individual.getGene1()[i] = temp3.get(i);
				} else {
					individual.getGene1()[i] = presentVertex[k];
					k++;
				}
			}
			population.add(individual);
		}
		return population;
	}

}
