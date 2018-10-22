package filesinout;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import structures.Cluster;
import structures.Vertex;

public class ReadFiles {
	public static ArrayList<Cluster> clusters = new ArrayList<Cluster>();
	public static Vertex[] vertices = new Vertex[5000];
	public static int num_vertex;
	// public static double weightMatrix[][] = readWeightMatrix();
	public static double weightMatrix[][];
	public static int numberOfCluster;
	public static int root;

	// instance 2
	public static ArrayList<Cluster> clusters1 = new ArrayList<Cluster>();
	public static Vertex[] vertices1 = new Vertex[5000];
	public static int num_vertex1;
	// public static double weightMatrix[][] = readWeightMatrix();
	public static double weightMatrix1[][];
	public static int numberOfCluster1;
	public static int root1;

	// public static int clusterReadFiles(){
	public static int clusterReadFiles(String fileName) {
		clusters.clear();
		BufferedReader br = null; // string to store data from file
		try {
			String sCurrentLine = null;
			br = new BufferedReader(new FileReader(fileName));
			// read lines 1..4
			for (int j = 0; j < 3; j++) {
				sCurrentLine = br.readLine();
			}
			// System.out.println(" str = "+str[1]);
			String[] str = sCurrentLine.split(": ");
			// System.out.println(" str = "+str[1]);
			num_vertex = Integer.parseInt(str[1]);

			// weightMatrix = new double[num_vertex][num_vertex];
			weightMatrix = new double[num_vertex][num_vertex];
			sCurrentLine = br.readLine();
			str = sCurrentLine.split(": ");
			numberOfCluster = Integer.parseInt(str[1]);
			sCurrentLine = br.readLine();
			sCurrentLine = br.readLine();

			// read the detail of the vertex
			for (int j = 0; j < num_vertex; j++) {
				sCurrentLine = br.readLine();
				str = sCurrentLine.split("\\s+");
				vertices[j] = new Vertex();
				// set coordinates to city
				vertices[j].setX(Double.parseDouble(str[1]));
				vertices[j].setY(Double.parseDouble(str[2]));
				// calculate distances
				for (int i = 0; i <= j; i++) {
					if (i == j) {
						weightMatrix[j][i] = 0;
					} else {
						weightMatrix[j][i] = weightMatrix[i][j] = Math
								.sqrt(Math.pow((vertices[j].getX() - vertices[i].getX()), 2)
										+ Math.pow((vertices[j].getY() - vertices[i].getY()), 2));
					}
				}
			}

			sCurrentLine = br.readLine();
			sCurrentLine = br.readLine();
			str = sCurrentLine.split(": ");
			root = Integer.parseInt(str[1]);
			for (int i = 0; i < numberOfCluster; i++) {
				int numberClusterVertex;
				int arrayCluster;
				sCurrentLine = br.readLine();
				str = sCurrentLine.split(" ");
				numberClusterVertex = str.length;
				Cluster cluster = new Cluster();
				for (int j = 0; j < numberClusterVertex - 2; ++j) {
					arrayCluster = Integer.parseInt(str[j + 1]);
					cluster.addElement(arrayCluster, j);
				}
				// System.out.println();
				clusters.add(cluster);
			}

		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		return num_vertex; // return number of vertex
	}

	public static int clusterReadFiles1(String fileName) {
		clusters1.clear();
		BufferedReader br = null; // string to store data from file
		try {
			String sCurrentLine = null;
			br = new BufferedReader(new FileReader(fileName));
			// read lines 1..4
			for (int j = 0; j < 3; j++) {
				sCurrentLine = br.readLine();
			}
			// System.out.println(" str = "+str[1]);
			String[] str = sCurrentLine.split(": ");
			// System.out.println(" str = "+str[1]);
			num_vertex1 = Integer.parseInt(str[1]);

			// weightMatrix = new double[num_vertex][num_vertex];
			weightMatrix1 = new double[num_vertex1][num_vertex1];
			sCurrentLine = br.readLine();
			str = sCurrentLine.split(": ");
			numberOfCluster1 = Integer.parseInt(str[1]);
			sCurrentLine = br.readLine();
			sCurrentLine = br.readLine();

			// read the detail of the vertex
			for (int j = 0; j < num_vertex1; j++) {
				sCurrentLine = br.readLine();
				str = sCurrentLine.split("\\s+");
				vertices1[j] = new Vertex();
				// set coordinates to city
				vertices1[j].setX(Double.parseDouble(str[1]));
				vertices1[j].setY(Double.parseDouble(str[2]));
				// calculate distances
				for (int i = 0; i <= j; i++) {
					if (i == j) {
						weightMatrix1[j][i] = 0;
					} else {
						weightMatrix1[j][i] = weightMatrix1[i][j] = Math
								.sqrt(Math.pow((vertices1[j].getX() - vertices1[i].getX()), 2)
										+ Math.pow((vertices1[j].getY() - vertices1[i].getY()), 2));
					}
				}
			}

			sCurrentLine = br.readLine();
			sCurrentLine = br.readLine();
			str = sCurrentLine.split(": ");
			root1 = Integer.parseInt(str[1]);
			for (int i = 0; i < numberOfCluster1; i++) {
				int numberClusterVertex;
				int arrayCluster;
				sCurrentLine = br.readLine();
				str = sCurrentLine.split(" ");
				numberClusterVertex = str.length;
				Cluster cluster = new Cluster();
				for (int j = 0; j < numberClusterVertex - 2; ++j) {
					arrayCluster = Integer.parseInt(str[j + 1]);
					cluster.addElement(arrayCluster, j);
				}
				// System.out.println();
				clusters1.add(cluster);
			}

		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		return num_vertex1; // return number of vertex
	}

	public static double[][] readWeightMatrix() {
		String filename = "C:/Users/TrungTB/Desktop/ClusteredTree_NON_EUC_Fom_Clustered_Tree_FINAL/Type_1_Small/5eil51.clt";
		Scanner sc;
		double[][] weghtMatrix2 = new double[num_vertex][num_vertex];
		try {
			sc = new Scanner(new File(filename));
			sc.nextLine();
			sc.nextLine();
			sc.nextLine();
			sc.nextLine();
			sc.nextLine();
			sc.nextLine();
			// optimalCost = sc.nextDouble();
			for (int i = 0; i < num_vertex; i++) {
				for (int j = 0; j < num_vertex; j++) {
					weghtMatrix2[i][j] = sc.nextDouble();
					System.out.print(" " + weghtMatrix2[i][j]);
				}
				System.out.println();
			}

		} catch (FileNotFoundException e) {

			e.printStackTrace();
		}

		return weghtMatrix2;

	}
}
