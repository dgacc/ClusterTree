package Cluster_Tree;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class Individual {
	public  static int num_vertex = ReadFiles.num_vertex;
	public  double[][] gene =  new double[num_vertex][num_vertex]; 
	InitializeChromosome initializechromosome = new InitializeChromosome();
	Evaluation evaluation = new Evaluation();
	private double[] constraintViolation;
	private double[] factorialCost;
	double cost;
	private int[] factorialRank = new int[2];
	double scalarFitness;
	int skillFactor;
	
	
	public Individual() {
	}
    public void inintilizeIndividual(){
    	gene = initializechromosome.clusterPrimRST(ReadFiles.weightMatrix, ReadFiles.clusters, num_vertex);
 
    }
    // getter
    public double[][] getGene(){
    	return gene;
    }
    //setter
    public void setGene(double[][] gene){
    	 this.gene = gene;	
    }
    
    // get fitness
    public double getFitness(Individual individual){
    	 return evaluation.evaluation(ReadFiles.weightMatrix, individual.getGene(), num_vertex, ReadFiles.root);
    }
	public void setFactorialRank(int rank, int index) {
		this.factorialRank[index] = rank;
	}

	public double[] getConstraintViolation() {
		return constraintViolation;
	}

	public void setConstraintViolation(double[] constraintViolation) {
		this.constraintViolation = constraintViolation;
	}

	public double[] getFactorialCost() {
		return factorialCost;
	}

	public void setFactorialCost(double[] factorialCost) {
		this.factorialCost = factorialCost;
	}

	public int[] getFactorialRank() {
		return factorialRank;
	}

	public void setFactorialRank(int[] factorialRank) {
		this.factorialRank = factorialRank;
	}

	public double getScalarFitness() {
		return scalarFitness;
	}

	public void setScalarFitness(double scalarFitness) {
		this.scalarFitness = scalarFitness;
	}

	public int getSkillFactor() {
		return skillFactor;
	}

	public void setSkillFactor(int skillFactor) {
		this.skillFactor = skillFactor;
	}
   
}