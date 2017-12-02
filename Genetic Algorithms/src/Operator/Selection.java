package Operator;

import java.util.Random;

import Structures.Individual;
import Structures.Population;

public class Selection {
	Individual individual  = new Individual();
	
	
	public int  touranmentSelection(double[] populationfitness, int tournament_size, Random r) {
		int[] tournamentGroup = new int[tournament_size];
		int index = 0;
		for(int i = 0 ; i < tournament_size; i++){
		tournamentGroup[i] = r.nextInt(Population.populationLength);	
		}
		
		double tempTour = populationfitness[tournamentGroup[0]];
		
		for( int i = 1; i < tournament_size; i++ ){
			if( tempTour < populationfitness[tournamentGroup[i]] ){
				tempTour = populationfitness[tournamentGroup[i]];
				index = tournamentGroup[i];
			}
		}
		return index;
		
	}


}

