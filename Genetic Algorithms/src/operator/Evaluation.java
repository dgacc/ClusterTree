package operator;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;
import filesinout.ReadFiles;
import structures.Cluster;
import structures.Population;

public class Evaluation {
	GraphMethods graphMethodsClass = new GraphMethods();
	CayleyCode cayleyCode = new CayleyCode();

	/**
	 * evaluate the fitness for cluster tree
	 * 
	 * @param weightMatrix
	 * @param edgeMatrix
	 * @param num_vertex
	 * @param startVertex
	 * @return pathLength: meaning the distance from root to the others vertex;
	 */
	public double clusterEvaluate(double[][] weightMatrix, double[][] edgeMatrix, int num_vertex, int startVertex) {
		double pathLength = 0;
		int[] pre = new int[num_vertex];
		ArrayList<Integer> path = new ArrayList<Integer>();
		double[][] tempMatrix = new double[num_vertex][num_vertex];

		// initialize the value to new weight matrix
		for (int i = 0; i < num_vertex; i++) {
			for (int j = 0; j < num_vertex; j++) {
				if (edgeMatrix[i][j] > 0) {
					tempMatrix[i][j] = weightMatrix[i][j];
				} else {
					tempMatrix[i][j] = 0;
				}
			}
		}
		// start to the distances;
		pre = graphMethodsClass.dijkstra(weightMatrix, num_vertex, startVertex);
		for (int i = 1; i < num_vertex; i++) {
			path = graphMethodsClass.printPath(startVertex, i, pre);
			for (int j = path.size() - 1; j > 0; j--) {
				pathLength = pathLength + tempMatrix[path.get(j)][path.get(j - 1)];

			}
		}
		return pathLength;
	}

	public double evaluation(double[][] weightMatrix, double[][] tree, int num_vertex, int startVertex) {
		double[] distances = new double[num_vertex];// distance between root and
													// the others
		double sum = 0;
		distances[startVertex] = 0;
		boolean[] mark = new boolean[num_vertex];
		Queue<Integer> queue = new LinkedList<>();
		for (int i = 0; i < num_vertex; i++) {
			mark[i] = true;
		}
		queue.add(startVertex);
		while (!queue.isEmpty()) {
			int u = queue.poll();
			mark[u] = false;
			for (int i = 0; i < num_vertex; i++) {
				if (tree[u][i] > 0 && mark[i]) {
					queue.add(i);
					mark[i] = false;
					distances[i] = distances[u] + weightMatrix[u][i];
					sum += distances[i];

				}
			}
		}
		return sum;
	}

	/**
	 * 
	 * @param weightMatrix
	 * @param tree
	 * @param num_vertex
	 * @param startVertex
	 * @return
	 */
	public double distanceEvaluate(double[][] weightMatrix, double[][] tree, int num_vertex, int startVertex) {
		int considerVertex;
		double sum = 0;
		boolean[] mark = new boolean[num_vertex];
		Queue<Integer> queue = new LinkedList<>();
		// initialize for : distances root, mark = true, and etc.....
		for (int i = 0; i < num_vertex; i++) {
			mark[i] = true;
		}
		queue.add(startVertex);
		mark[startVertex] = false;
		while (!queue.isEmpty()) {
			considerVertex = queue.poll();
			for (int i = 0; i < num_vertex; i++) {
				if (tree[considerVertex][i] > 0 && mark[i]) {
					sum = sum + weightMatrix[considerVertex][i];
					queue.add(i);
					mark[i] = false;
				}
			}

		}
		return sum;
	}

	/*********************************************************************************************************************************************
	 * Decoding cho bai toan cay khung Giam tu n dinh ve thanh m dinh 01. Xóa
	 * các đỉnh và cạnh liên thuộc với nó > m 02. Tìm các thành phần liên thông
	 * 03. Nối các thành phần liên thông thứ i -> i + 1 + Ch�?n ngẫu nhiên ở mỗi
	 * thành phần liên thông 1 đỉnh + Ch�?n đỉnh đầu tiên
	 ********************************************************************************************************************************************/
	public double[][] decodingMFOVertexInSubGraph(double[][] ind_Matrix, int max_Genes, int num_Gen_of_Task_j) {
		double[][] tmp_Matrix = new double[max_Genes][max_Genes];
		// 01. Tìm các đỉnh lớn hơn num_Gen_of_Task_j và xóa nó cùng cạnh liên
		// thuộc
		for (int i = 0; i < num_Gen_of_Task_j; i++) {
			for (int j = 0; j < num_Gen_of_Task_j; j++) {
				if ((ind_Matrix[i][j] > 0) && (ind_Matrix[i][j] < Double.MAX_VALUE)) {
					tmp_Matrix[i][j] = 1.0f;
				} else {
					tmp_Matrix[i][j] = 0.0f;
				}
			}
		}
		// 02. Tìm các thành phần liên thông va lay moi tp lien thong 1 dinh
		int[] tp_LT = graphMethodsClass.get_Vertex_In_Each_SubGraph(tmp_Matrix, num_Gen_of_Task_j);

		// 03. Nối các thành phần liên thông thứ i -> i + 1
		for (int i = 0; i < tp_LT.length - 1; i++) {
			tmp_Matrix[tp_LT[i]][tp_LT[i + 1]] = 1.0f;
			tmp_Matrix[tp_LT[i + 1]][tp_LT[i]] = 1.0f;
		}

		// Tao ra ma tran ket qua
		double[][] final_Matrix = new double[num_Gen_of_Task_j][num_Gen_of_Task_j];
		for (int i = 0; i < num_Gen_of_Task_j; i++) {
			for (int j = 0; j < num_Gen_of_Task_j; j++) {
				final_Matrix[i][j] = (int) tmp_Matrix[i][j];
			}
		}
		return final_Matrix;
	}

	/**
	 * In order to make tree by using the prufer number code decoding.
	 * 
	 * @param pruferNumberCode
	 *            : the array of prufer number code
	 * @param num_vertices
	 *            : the number of vertices
	 * @return tree after apply make tree prufer number algorithm to prufer
	 *         number code
	 */
	public double[][] pruferNumber(int[] pruferNumberCode, int num_vertices) {

		double[][] tree = new double[num_vertices][num_vertices];
		for (int i = 0; i < num_vertices; i++) {
			for (int j = 0; j < num_vertices; j++)
				tree[i][j] = -0.001;
		}
		int[] degree = new int[num_vertices]; // store the degree of vertices
		int[] tempTable = new int[num_vertices]; // store the node from 1 to n

		for (int i = 0; i < num_vertices; i++) {
			degree[i] = 1; // initialize the degree of vertices = 1
			tempTable[i] = i; // temptable = vertices from 1 to n
		}

		for (int i = 0; i < num_vertices - 2; i++) {
			int temp = pruferNumberCode[i];
			degree[temp]++;
		}

		// build the tree from prufer number code
		for (int i = 0; i < num_vertices - 2; i++) {
			int temp = pruferNumberCode[i];
			for (int j = 0; j < num_vertices; j++) {
				if (degree[j] == 1) {
					tree[temp][j] = tree[j][temp] = 1.0;
					degree[temp]--;
					degree[j]--;
					break;
				}
			}
		}
		// if prufer number code have no element left
		int flag1 = 0;
		int flag2 = 0;
		for (int i = 0; i < num_vertices; i++) {

			if ((degree[i] == 1) && (flag1 == 0)) {
				flag1 = i;
			} else if ((degree[i] == 1)) {
				flag2 = i;
				break;
			}

		}
		tree[flag1][flag2] = tree[flag2][flag1] = 1.0f;
		return tree;
	}

	/**
	 * 
	 * @param pruferNumberCode
	 * @param startpoint
	 * @param clusters
	 * @param maxClusters
	 * @param num_vertex
	 * @param rnd
	 * @return
	 */

	public double[][] decodingPruferNumber(int[] pruferNumberCode, int[] startpoint, ArrayList<Cluster> clusters,
			ArrayList<Cluster> maxClusters, int num_vertex, Random rnd) {
		int numberOfCluster = clusters.size();
		int maxNumberOfCluster = maxClusters.size();

		double[][] tree = new double[num_vertex][num_vertex];
		int[] indexCluster = new int[numberOfCluster];
		int count1 = 0;
		int[] presentationClusterPrufer = new int[numberOfCluster - 2];
		for (int i = 0; i < numberOfCluster; i++) {
			// lay cho cluster dai dien

			if (i < numberOfCluster - 2 && pruferNumberCode[startpoint[maxNumberOfCluster - 1] + i] < numberOfCluster) {
				// presentationClusterPrufer[i] =
				// pruferNumberCode[startpoint[maxNumberOfCluster - 1] + i];
				presentationClusterPrufer[count1] = pruferNumberCode[startpoint[maxNumberOfCluster - 1] + i];
				count1++;
			}

			int numberClusterVertices = clusters.get(i).getCluster().size();
			double[][] clusterTree = new double[numberClusterVertices][numberClusterVertices];

			int[] pruferNumberCodeCluster = new int[numberClusterVertices - 2];
			// get prufer number code of each cluster
			int count = 0;
			for (int j = 0; j < numberClusterVertices - 2; j++) {
				if (i == 0) {

					if (pruferNumberCode[j] < numberClusterVertices) {
						pruferNumberCodeCluster[count] = pruferNumberCode[j];
						count++;
					}

				} else {
					if (pruferNumberCode[startpoint[i - 1] + j] < numberClusterVertices) {
						pruferNumberCodeCluster[count] = pruferNumberCode[startpoint[i - 1] + j];
						count++;
					}
				}
			}
			int temp1 = (numberClusterVertices - count - 2);
			for (int k = 0; k < temp1; k++) {
				int temp = rnd.nextInt(numberClusterVertices);
				pruferNumberCodeCluster[k + count] = temp;
			}
			// build tree from prufer number code
			clusterTree = pruferNumber(pruferNumberCodeCluster, numberClusterVertices);
			// clusterTree =
			// cayleyCode.decodingDandelion(pruferNumberCodeCluster);
			// clusterTree = cayleyCode.decodingBlob(pruferNumberCodeCluster);
			for (int j = 0; j < numberClusterVertices; j++) {
				for (int k = 0; k < numberClusterVertices; k++) {
					tree[clusters.get(i).getCluster().get(j)][clusters.get(i).getCluster().get(k)] = clusterTree[j][k];
				}
			}
			int position = pruferNumberCode[startpoint[maxNumberOfCluster - 1] + maxNumberOfCluster - 2 + i];
			if (position >= numberClusterVertices) {
				position = position % numberClusterVertices;
			}
			indexCluster[i] = clusters.get(i).getCluster().get(position);

		}
		int temp2 = (numberOfCluster - count1 - 2);
		for (int k = 0; k < temp2; k++) {
			int temp = rnd.nextInt(numberOfCluster);
			presentationClusterPrufer[k + count1] = temp;
		}

		double[][] spanningTreeOfCluster = pruferNumber(presentationClusterPrufer, numberOfCluster);
		// double[][] spanningTreeOfCluster =
		// cayleyCode.decodingDandelion(presentationClusterPrufer);
		// double[][] spanningTreeOfCluster =
		// cayleyCode.decodingBlob(presentationClusterPrufer);
		for (int i = 0; i < numberOfCluster; i++) {
			for (int j = 0; j < numberOfCluster; j++) {
				tree[indexCluster[i]][indexCluster[j]] = spanningTreeOfCluster[i][j];
			}
		}
		return tree;

	}

	/**
	 * In order to build up information for dividing prufer number code to
	 * prufer number code for each cluster.
	 * 
	 * @param maxClusters
	 * @return array of start point: each element is value of start point of
	 *         prufer number of each cluster
	 */

	public int[] getStartPoint(ArrayList<Cluster> maxClusters) {
		int numberOfMaxCluster = maxClusters.size();
		int[] startpoint = new int[numberOfMaxCluster];// start point to divide
														// the prufer number
														// code to build tree
														// for each cluster.
		if (maxClusters.get(0).getCluster().size() < 3) {
			startpoint[0] = 1;
		} else {
			startpoint[0] = maxClusters.get(0).getCluster().size() - 2;
		}
		for (int i = 1; i < numberOfMaxCluster; i++) {
			if (maxClusters.get(i).getCluster().size() < 3) {
				startpoint[i] = startpoint[i - 1] + 1;
			} else {
				startpoint[i] = startpoint[i - 1] + maxClusters.get(i).getCluster().size() - 2;
			}
		}
		return startpoint;
	}

	/**
	 * 
	 * @param pruferNumberCode
	 * @param startpoint
	 * @param clusters
	 * @param num_vertex
	 * @return
	 */
	public double[][] decodingPruferNumber(int[] pruferNumberCode, int[] startpoint, ArrayList<Cluster> clusters,
			int num_vertex) {

		int numberOfCluster = clusters.size();
		double[][] tree = new double[num_vertex][num_vertex];
		int[] indexCluster = new int[numberOfCluster];
		int[] presentationClusterPrufer = new int[numberOfCluster - 2];
		for (int i = 0; i < numberOfCluster; i++) {
			if (i < numberOfCluster - 2) {
				presentationClusterPrufer[i] = pruferNumberCode[startpoint[numberOfCluster - 1] + i];
			}
			int numberClusterVertices = clusters.get(i).getCluster().size();
			double[][] clusterTree = new double[numberClusterVertices][numberClusterVertices];

			int[] pruferNumberCodeCluster = new int[numberClusterVertices - 2];
			// get prufer number code of each cluster
			for (int j = 0; j < numberClusterVertices - 2; j++) {
				if (i == 0) {
					pruferNumberCodeCluster[j] = pruferNumberCode[j];
				} else {
					pruferNumberCodeCluster[j] = pruferNumberCode[startpoint[i - 1] + j];

				}
			}
			// build tree from prufer number code
			// clusterTree = pruferNumber(pruferNumberCodeCluster,
			// numberClusterVertices);
			clusterTree = cayleyCode.decodingDandelion(pruferNumberCodeCluster);
			for (int j = 0; j < numberClusterVertices; j++) {
				for (int k = 0; k < numberClusterVertices; k++) {
					tree[clusters.get(i).getCluster().get(j)][clusters.get(i).getCluster().get(k)] = clusterTree[j][k];
				}
			}
			int position = pruferNumberCode[startpoint[numberOfCluster - 1] + numberOfCluster - 2 + i];
			indexCluster[i] = clusters.get(i).getCluster().get(position);

		}

		// double[][] spanningTreeOfCluster =
		// pruferNumber(presentationClusterPrufer, numberOfCluster);
		double[][] spanningTreeOfCluster = cayleyCode.decodingDandelion(presentationClusterPrufer);

		for (int i = 0; i < numberOfCluster; i++) {
			for (int j = 0; j < numberOfCluster; j++) {
				tree[indexCluster[i]][indexCluster[j]] = spanningTreeOfCluster[i][j];
			}
		}
		// Windows windows = new Windows();
		// windows.runWindow(" cac cay");
		// Paint p = new Paint();
		// p.setPaint(tree, ReadFiles.vertices1, ReadFiles.clusters1,
		// ReadFiles.num_vertex1,
		// 0, 0, ReadFiles.root1);
		// windows.addPaint(p);
		return tree;
	}

	public double[] getPopFitness(Population population, int populationLength) {
		double[] popFitness = new double[populationLength];
		for (int i = 0; i < populationLength; i++) {
			double[][] tree = decodingPruferNumber(population.getIndividual(i).getGene1(),
					getStartPoint(ReadFiles.clusters), ReadFiles.clusters, ReadFiles.num_vertex);

			popFitness[i] = evaluation(ReadFiles.weightMatrix, tree, ReadFiles.num_vertex, ReadFiles.root);
		}
		return popFitness;
	}

	public double[][] decodingCayley(int[] pruferNumberCode, int[] startpoint, ArrayList<Cluster> clusters,
			ArrayList<Cluster> maxClusters, int num_vertex, Random rnd) {
		int numberOfCluster = clusters.size();
		int maxNumberOfCluster = maxClusters.size();

		double[][] tree = new double[num_vertex][num_vertex];
		int[] indexCluster = new int[numberOfCluster];
		int count1 = 0;
		int count4 = 0;
		int[] presentationClusterPrufer = new int[numberOfCluster - 2];
		int[] auxilityArray1 = new int[numberOfCluster - 2];
		for (int i = 0; i < numberOfCluster; i++) {
			// lay cho cluster dai dien

			if (i < numberOfCluster - 2) {
				if (pruferNumberCode[startpoint[maxNumberOfCluster - 1] + i] < numberOfCluster) {
					presentationClusterPrufer[count1] = pruferNumberCode[startpoint[maxNumberOfCluster - 1] + i];
					count1++;
				}
			} else {
				auxilityArray1[count4] = pruferNumberCode[startpoint[maxNumberOfCluster - 1] + i] % numberOfCluster;
				count4++;
			}
			int numberClusterVertices = clusters.get(i).getCluster().size();
			double[][] clusterTree = new double[numberClusterVertices][numberClusterVertices];

			int[] pruferNumberCodeCluster = new int[numberClusterVertices - 2];
			int[] auxilityArray = new int[numberClusterVertices - 2];
			// get prufer number code of each cluster
			int count = 0;
			int count2 = 0;
			for (int j = 0; j < numberClusterVertices - 2; j++) {
				if (i == 0) {

					if (pruferNumberCode[j] < numberClusterVertices) {
						pruferNumberCodeCluster[count] = pruferNumberCode[j];
						count++;
					} else {
						auxilityArray[count2] = pruferNumberCode[j] % numberClusterVertices;
						count2++;
					}

				} else {
					if (pruferNumberCode[startpoint[i - 1] + j] < numberClusterVertices) {
						pruferNumberCodeCluster[count] = pruferNumberCode[startpoint[i - 1] + j];
						count++;
					} else {
						auxilityArray[count2] = pruferNumberCode[startpoint[i - 1] + j] % numberClusterVertices;
						count2++;
					}
				}
			}
			int temp1 = (numberClusterVertices - count - 2);
			for (int k = 0; k < temp1; k++) {
				pruferNumberCodeCluster[k + count] = auxilityArray[k];
			}
			// build tree from prufer number code
			// clusterTree = pruferNumber(pruferNumberCodeCluster,
			// numberClusterVertices);
			clusterTree = cayleyCode.decodingDandelion(pruferNumberCodeCluster);
			// clusterTree = cayleyCode.decodingBlob(pruferNumberCodeCluster);
			for (int j = 0; j < numberClusterVertices; j++) {
				for (int k = 0; k < numberClusterVertices; k++) {
					tree[clusters.get(i).getCluster().get(j)][clusters.get(i).getCluster().get(k)] = clusterTree[j][k];
				}
			}
			int position = pruferNumberCode[startpoint[maxNumberOfCluster - 1] + maxNumberOfCluster - 2 + i];
			if (position >= numberClusterVertices) {
				position = position % numberClusterVertices;
			}
			indexCluster[i] = clusters.get(i).getCluster().get(position);

		}
		int temp2 = (numberOfCluster - count1 - 2);
		for (int k = 0; k < temp2; k++) {
			int temp = rnd.nextInt(numberOfCluster);
			presentationClusterPrufer[k + count1] = temp;
		}

		// double[][] spanningTreeOfCluster =
		// pruferNumber(presentationClusterPrufer, numberOfCluster);
		double[][] spanningTreeOfCluster = cayleyCode.decodingDandelion(presentationClusterPrufer);
		// double[][] spanningTreeOfCluster =
		// cayleyCode.decodingBlob(presentationClusterPrufer);
		for (int i = 0; i < numberOfCluster; i++) {
			for (int j = 0; j < numberOfCluster; j++) {
				tree[indexCluster[i]][indexCluster[j]] = spanningTreeOfCluster[i][j];
			}
		}
		return tree;

	}

	/**
	 * In order to decode Chromosome in unified search space into solution for a
	 * task
	 * 
	 * @param Chromosome
	 *            Chromosome in unified search space
	 * @param startpoint
	 *            The point separate all of segments
	 * @param clusters
	 *            The list of clusters of the problem which is decoded from
	 *            chromosome unified search space
	 * @param maxClusters
	 *            The list of cluster of Individual in unified search space
	 * @param num_vertex
	 *            Number of vertices in Problem which is decoded from chromosome
	 *            unified search space
	 * @param rnd
	 *            Random state
	 * @return Tree(solution) after apply decoding method.
	 */
	public double[][] decodingCayleyTest(int[] Chromosome, int[] startpoint, ArrayList<Cluster> clusters,
			ArrayList<Cluster> maxClusters, int num_vertex, Random rnd) {
		int numberOfCluster = clusters.size();
		int maxNumberOfCluster = maxClusters.size();

		double[][] tree = new double[num_vertex][num_vertex];
		int[] indexCluster = new int[numberOfCluster];
		int count1 = 0;
		int count4 = 0;
		int[] presentationClusterPrufer = new int[numberOfCluster - 2];
		int[] auxilityArray1 = new int[numberOfCluster - 2];
		for (int i = 0; i < numberOfCluster; i++) {
			// lay cho cluster dai dien

			if (i < numberOfCluster - 2) {
				if (Chromosome[startpoint[maxNumberOfCluster - 1] + i] < numberOfCluster) {
					presentationClusterPrufer[count1] = Chromosome[startpoint[maxNumberOfCluster - 1] + i];
					count1++;
				}
			} else {
				auxilityArray1[count4] = Chromosome[startpoint[maxNumberOfCluster - 1] + i] % numberOfCluster;
				count4++;
			}
			int numberClusterVertices = clusters.get(i).getCluster().size();
			double[][] clusterTree = new double[numberClusterVertices][numberClusterVertices];
			// check whether number elements of each cluster greater 2 or not?
			if (numberClusterVertices == 1) {
				// do nothing
			} else if (numberClusterVertices == 2) {
				clusterTree[0][1] = 1.0f;
				clusterTree[1][0] = 1.0f;

			} else {

				int[] pruferNumberCodeCluster = new int[numberClusterVertices - 2];
				int[] auxilityArray = new int[numberClusterVertices - 2];
				// get Cayley code for each cluster
				int count = 0;
				int count2 = 0;
				for (int j = 0; j < numberClusterVertices - 2; j++) {
					if (i == 0) {

						if (Chromosome[j] < numberClusterVertices) {
							pruferNumberCodeCluster[count] = Chromosome[j];
							count++;
						} else {
							auxilityArray[count2] = Chromosome[j] % numberClusterVertices;
							count2++;
						}

					} else {
						if (Chromosome[startpoint[i - 1] + j] < numberClusterVertices) {
							pruferNumberCodeCluster[count] = Chromosome[startpoint[i - 1] + j];
							count++;
						} else {
							auxilityArray[count2] = Chromosome[startpoint[i - 1] + j] % numberClusterVertices;
							count2++;
						}
					}
				}
				int temp1 = (numberClusterVertices - count - 2);
				for (int k = 0; k < temp1; k++) {
					pruferNumberCodeCluster[k + count] = auxilityArray[k];
				}
				// build tree from prufer number code
				// clusterTree = pruferNumber(pruferNumberCodeCluster,
				// numberClusterVertices);
//				clusterTree = cayleyCode.decodingDandelion(pruferNumberCodeCluster);
				 clusterTree = cayleyCode.decodingBlob(pruferNumberCodeCluster);
			}
			for (int j = 0; j < numberClusterVertices; j++) {
				for (int k = 0; k < numberClusterVertices; k++) {
					tree[clusters.get(i).getCluster().get(j)][clusters.get(i).getCluster().get(k)] = clusterTree[j][k];
				}
			}
			int position = Chromosome[startpoint[maxNumberOfCluster - 1] + maxNumberOfCluster - 2 + i];
			if (position >= numberClusterVertices) {
				position = position % numberClusterVertices;
			}
			indexCluster[i] = clusters.get(i).getCluster().get(position);

		}
		int temp2 = (numberOfCluster - count1 - 2);
		for (int k = 0; k < temp2; k++) {
			int temp = rnd.nextInt(numberOfCluster);
			presentationClusterPrufer[k + count1] = temp;
		}

		// double[][] spanningTreeOfCluster =
		// pruferNumber(presentationClusterPrufer, numberOfCluster);
//		double[][] spanningTreeOfCluster = cayleyCode.decodingDandelion(presentationClusterPrufer);
		 double[][] spanningTreeOfCluster = cayleyCode.decodingBlob(presentationClusterPrufer);
		for (int i = 0; i < numberOfCluster; i++) {
			for (int j = 0; j < numberOfCluster; j++) {
				tree[indexCluster[i]][indexCluster[j]] = spanningTreeOfCluster[i][j];
			}
		}
		return tree;

	}
}
