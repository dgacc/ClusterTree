package Dislay;

import javax.swing.JFrame;

public class Windows {
	public static Paint  p = new Paint();
	
	public static void  main(String[] ags){
		JFrame gf = new JFrame();
		gf.setVisible(true);
		gf.setSize(480, 720);
		gf.setTitle("TrungTB");
		gf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		// add new painting in here
	 
		gf.add(p);
		gf.setVisible(true);
		}
}
