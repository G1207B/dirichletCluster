package edu.ucr.cuilab.algorithm.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import edu.ucr.cuilab.algorithm.EncodeBinTask;

public class BuildTestData {

	public static void overlapSample() throws Exception {
		String originFile = EncodeBinTask.FILENAME;
		String overlapInfo = EncodeBinTask.OUTPUT;
		String outputFile = originFile + ".sample";
		List<String> strList = new ArrayList<String>();
		Set<Integer> selectSet = new TreeSet<Integer>();
		BufferedReader overlapBr = new BufferedReader(new FileReader(new File(overlapInfo)));
		String line = null;
		while (null != (line = overlapBr.readLine())) {
			String[] split = line.split("\\s\\s*");
			for (String info:split) {
				try {
					selectSet.add(Integer.valueOf(info));
				} catch (NumberFormatException nfe) {
					continue;
				}
			}
		}
		overlapBr.close();
		
		BufferedReader originBr = new BufferedReader(new FileReader(new File(originFile)));
		while (null != (line = originBr.readLine())) {
			strList.add(line);
		}
		originBr.close();
		
		PrintWriter pw = new PrintWriter(new FileWriter(new File(outputFile)));
		for (Integer num:selectSet) {
			pw.println(strList.get(num * 2));
			pw.println(strList.get(num * 2 + 1));
		}
		pw.close();
		
	}
	
	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub

		overlapSample();
	}

}
