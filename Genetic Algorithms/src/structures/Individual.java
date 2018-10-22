package structures;

import filesinout.ReadFiles;
import operator.Evaluation;
import operator.InitializeChromosome;

public class Individual {
	public static int num_vertex = ReadFiles.num_vertex;
	private double[][] gene = new double[num_vertex][num_vertex];
	private int[] gene1 = new int[5000];
	InitializeChromosome initializechromosome = new InitializeChromosome();
	Evaluation evaluation = new Evaluation();
	protected double[] constraintViolation;
	protected double[] factorialCost;
	public double cost;
	protected int[] factorialRank = new int[2];
	protected double scalarFitness;
	protected int skillFactor;

	public Individual() {
	}

	public void inintilizeIndividual() {
		this.gene = initializechromosome.clusterPrimRST(ReadFiles.weightMatrix, ReadFiles.clusters,
				ReadFiles.num_vertex);

	}

	// getter
	public double[][] getGene() {
		return this.gene;
	}

	// setter
	public void setGene(double[][] gene) {
		this.gene = gene;
	}

	// get fitness
	public double getFitness(Individual individual) {
		return evaluation.evaluation(ReadFiles.weightMatrix, individual.getGene(), ReadFiles.num_vertex,
				ReadFiles.root);
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

	public void setGene(int[] gene1) {
		this.gene1 = gene1;
	}

	public int[] getGene1() {
		return this.gene1;
	}

}