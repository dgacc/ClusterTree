package structures;

import operator.Evaluation;
import operator.InitializeChromosome;

public class Chromosome{
 
		InitializeChromosome initializechromosome = new InitializeChromosome();
		Evaluation evaluation = new Evaluation();
		protected double[] constraintViolation;
		protected double[] factorialCost;
		public double cost;
		protected int[] factorialRank = new int[2];
		protected double scalarFitness;
		protected int skillFactor;
		int[] gene;
		
		
	    // getter
	    public int[] getGene(){
	    	return this.gene;
	    }
	    //setter
	    public void setGene(int[] gene){
	    	 this.gene = gene;	
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