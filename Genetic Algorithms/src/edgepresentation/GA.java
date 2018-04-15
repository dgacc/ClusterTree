package edgepresentation;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import javax.swing.*;

import dislay.Paint;
import dislay.Windows;
import filesinout.ReadFiles;
import operator.Crossover;
import operator.Mutations;
import operator.Selection;
import random.MyRandom;
import structures.Individual;
import structures.Population;
public class GA{
	private static double m_rate = 0.05;
	private static final double p_rate = 0.5; // crossover rate
	private static  int generation = 500;
	public static int populationLength = 100; // default length of population
	private static Crossover crossover = new Crossover();
	private static Mutations mutation = new Mutations();
	private static Selection selection = new Selection();
/*

 */


	public static void main(String[] args) {
		String[] test = {"50pr439","75lin105","5i500-304","6i500","7i70-21","9pr439-3x3","2lin105-2x1","7i70-21"};
		for(int k= 0; k < test.length; k++){
		ReadFiles.clusterReadFiles("test/"+test[k]+".clt");
//	    ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/" +args[0]+ ".clt");
		int  num_vertex = ReadFiles.num_vertex;
//		 m_rate = 1.0/num_vertex;


//		Windows windows = new Windows();
//		windows.runWindow("The best ever");
		 String dirname = "Result/"+test[k];
	      File d1 = new File(dirname);
	      // Bay gio tao thu muc.
	      d1.mkdirs();

		 
//
	      System.out.println("-------------------------------------------------------------------------------");

		  for(int loop = 0; loop < 30; loop ++){
			MyRandom.setSeed(loop);
		    // write to files
		    PrintWriter pw = null;
			PrintWriter pw3 = null;
			
			
				try {
					pw = new PrintWriter(new FileWriter(dirname+"/"+test[k]+"-seed("+loop+").opt", true));
					pw3 = new PrintWriter(new FileWriter(dirname+"/"+test[k]+"_seed("+loop+").gen", true));
				} catch (IOException e) {
					e.printStackTrace();
				}
				pw3.println("Generations \t"+test[k]);
				
				 long start = System.currentTimeMillis();
			// Initialize population
		Population population = new Population();
		for( int i = 0; i < populationLength; i++  ){
	        Individual individual = new Individual();
	        individual.inintilizeIndividual();
			population.addIndiv(individual);
		}

		for (int i = 0; i <  generation; i++) {
			Population subPop = new Population(); // store new individuals are generated
			double[] popFitness = population.populationFitness();


			for (int j = 0; j < Population.populationLength / 2; j++) {
				Individual father = new  Individual();
				Individual mother = new  Individual();

				//use tournament to  select Individual
				father = population.getIndividual(selection.touranmentSelection(popFitness, 2, MyRandom.r));
				mother = population.getIndividual(selection.touranmentSelection(popFitness, 2, MyRandom.r));


				Individual offspring1 =  new Individual();
				Individual offspring2 =  new Individual();

				// new  random number  [0;1]
				double d = MyRandom.r.nextDouble();
				// if random number >  crossover rate,  add  the individual to population without crossover
		     	if (d < p_rate) {

			    	offspring1.setGene(crossover.clusterCrossover(father.getGene(), mother.getGene(),
	   		    		       num_vertex, ReadFiles.clusters, MyRandom.r));
			    	offspring2.setGene(crossover.clusterCrossover(father.getGene(), mother.getGene(),
		    		       num_vertex, ReadFiles.clusters,MyRandom.r));

		     	}else{
						 offspring1.setGene(mutation.mutationClusterTreeGA(mother.getGene(),
								 num_vertex, ReadFiles.clusters, MyRandom.r, m_rate));
						 offspring2.setGene(mutation.mutationClusterTreeGA(father.getGene(),
								 num_vertex, ReadFiles.clusters, MyRandom.r, m_rate));
					}


				subPop.addIndiv(offspring1);
				subPop.addIndiv(offspring2);

			}
		     	
		     	
		     	
			// get the index of the best Individual  of old Population
			double bestFiness = popFitness[0];
			int bestFinessIndex = 0;

			for(int t = 1; t < populationLength; t++ ){
			 if(popFitness[t] < bestFiness){
				 bestFiness = popFitness[t];
				  bestFinessIndex = t;
			 }

			}
			  pw3.println(i+"\t"+popFitness[bestFinessIndex]);
			// paint graph
//			Paint  p = new Paint();
//			p.setPaint(population.getIndividual(bestFinessIndex).getGene(), ReadFiles.vertices,
//					ReadFiles.clusters, num_vertex, popFitness[bestFinessIndex], 600, ReadFiles.root);
//			windows.addPaint(p);


			// add the best Individual  of old Population to new Population
			subPop.population.set(MyRandom.r.nextInt(Population.populationLength),
					population.getIndividual(bestFinessIndex));
					population = subPop;
		}
		pw3.close();
		long end = System.currentTimeMillis();
		double[] popFitness = population.populationFitness();
		double bestFiness = popFitness[0];
		int bestFinessIndex = 0;
		for(int t = 1; t < populationLength; t++ ){
			 if(popFitness[t] < bestFiness){
				 bestFiness = popFitness[t];
				  bestFinessIndex = t;
			 }
			}

	    // print to files
//		System.out.println(popFitness[bestFinessIndex]);
//		 pw.println(popFitness[bestFinessIndex]);
		MEF2Instances mef = new MEF2Instances();
		 System.out.println("|"+loop+"\t|"+test[k]+"\t|"+popFitness[bestFinessIndex]+"\t|"+mef.getDateFromMillis(end - start)+"\t|");
		 pw.println("Name: "+test[k]);
		 pw.println("Seed: " + loop);
		 pw.println("Fitness: " +popFitness[bestFinessIndex]);
		 pw.println("Time: " + mef.getDateFromMillis(end - start));
		 for(int m = 0; m < ReadFiles.num_vertex; m ++){
			 pw.println();
			 for( int n = 0; n < ReadFiles.num_vertex; n++){
				 if( population.getIndividual(bestFinessIndex).getGene()[m][n] > 0){
					 population.getIndividual(bestFinessIndex).getGene()[m][n] = 1.0f; }
				 pw.print((int) population.getIndividual(bestFinessIndex).getGene()[m][n]+"\t");
				
			 }
		 }
		 pw.close();

	}
		  }
}
}

