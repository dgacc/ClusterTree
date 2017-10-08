package Cluster_Tree;
import java.awt.*;
import java.util.ArrayList;
import java.util.Random;

import javax.swing.*;

class Paint extends JPanel {
    public int num_vertex;
    private static Vertex[] vertexs = ReadFiles.vertices;
    public double[][] weightMatrix ;
    private static ArrayList<Cluster> clusters = ReadFiles.clusters;
    public double fitness;
    
    
    
    
    private  void doDrawing(Graphics g) {
    	int numberOfCluster  = clusters.size();
        Graphics2D g2d = (Graphics2D) g;
        for( int i = 0; i < numberOfCluster; i ++ ){
        	int numberClusterVertex = clusters.get(i).getCluster().size();
        	Random r = new Random();
        	g2d.setColor(new Color(r.nextInt(256),r.nextInt(256),r.nextInt(256)));
        	for( int  j = 0; j < numberClusterVertex; j++ ){
        	g2d.fillOval((int)vertexs[clusters.get(i).getCluster().get(j)].getX()*10 - 5,
        			(int)vertexs[clusters.get(i).getCluster().get(j)].getY()*10 - 5, 10, 10);
        	}
        }
        for(int i = 0; i < num_vertex; i++){
        	for( int j = 0; j < num_vertex; j++){
        		 g2d.setColor(Color.RED);
        		 g2d.drawString(""+( i+ 1), (int)vertexs[i].getX()*10 + 10, (int)vertexs[i].getY()*10 + 10);
        		 if(weightMatrix[i][j] > 0){
        		 g2d.setColor(Color.BLUE);	 
                 g2d.drawLine((int)vertexs[i].getX()*10, (int)vertexs[i].getY()*10, 
                		 (int)vertexs[j].getX()*10, (int)vertexs[j].getY()*10);        		
        		 }
        }
   		}
        
        
        g2d.setColor(Color.GREEN);
        g2d.fillOval((int)vertexs[ReadFiles.root].getX()*10 - 10 ,
        		(int)vertexs[ReadFiles.root].getY()*10 - 10, 20, 20);
        g2d.drawString("fitness:" + fitness , 50, 20);
    }
 
    @Override
    public void paintComponent(Graphics g) {
 
        super.paintComponent(g);
        doDrawing(g);
    }
}
 
