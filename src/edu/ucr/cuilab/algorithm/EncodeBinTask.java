package edu.ucr.cuilab.algorithm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

public class EncodeBinTask {

	public static String FILENAME = "/home/xinping/Desktop/6008/abundance_species_equal.txt";
	public static String OUTPUT = "/home/xinping/Desktop/6008/encodeBin.test";

	public static void printUsage() {
		System.out
				.println("java EncodeBinTask abundance_species_equal.txt > output.bin");
		System.out.println("The sequences should be pair end");
	}

	public static String getTransReverse(String origin) {
		char[] charList = origin.toCharArray();
		for (int i = 0; i < charList.length; i++) {
			switch (charList[i]) {
			case 'A':
				charList[i] = 'T';
				break;
			case 'C':
				charList[i] = 'G';
				break;
			case 'G':
				charList[i] = 'C';
				break;
			case 'T':
				charList[i] = 'A';
				break;
			default:
				break;
			}
		}
		return (new StringBuilder(String.valueOf(charList))).reverse()
				.toString();
	}

	public static Set<Long> stringEncode(String seq, int wLen) {
		Set<Long> encodeSet = new HashSet<Long>();
		long encode0 = 0l;
		long encode1 = 0l;
		String tr = getTransReverse(seq);
		char[] seqCharArray = seq.toCharArray();
		char[] trCharArray = tr.toCharArray();
		int length = seq.length();

		for (int i = 0; i < wLen; i++) {
			encode0 <<= 2;
			encode1 <<= 2;
			switch (seqCharArray[i]) {
			case 'C':
				encode0 += 1;
				break;
			case 'G':
				encode0 += 2;
				break;
			case 'T':
				encode0 += 3;
			default:
			}

			switch (trCharArray[i]) {
			case 'C':
				encode1 += 1;
				break;
			case 'G':
				encode1 += 2;
				break;
			case 'T':
				encode1 += 3;
			default:
			}
		}

		encodeSet.add(encode0);
		encodeSet.add(encode1);

		for (int i = wLen; i < length; i++) {
			encode0 <<= 2;
			encode1 <<= 2;
			switch (seqCharArray[i]) {
			case 'C':
				encode0 += 1;
				break;
			case 'G':
				encode0 += 2;
				break;
			case 'T':
				encode0 += 3;
			default:
			}

			switch (trCharArray[i]) {
			case 'C':
				encode1 += 1;
				break;
			case 'G':
				encode1 += 2;
				break;
			case 'T':
				encode1 += 3;
			default:
			}
			encodeSet.add(encode0);
			encodeSet.add(encode1);
		}
		return encodeSet;
	}

	public static void mainJobTest(List<String> strList, int wLen)
			throws IOException {
		PrintWriter pw = new PrintWriter(new FileWriter(new File(OUTPUT)));

		for (int i = 0; i < strList.size(); i++) {
			Set<Long> encodeSet = stringEncode(strList.get(i), wLen);
			for (Long l : encodeSet) {
				pw.print("\t");
				pw.print(l);
			}
			pw.println();
		}
		pw.close();
	}
	
	

	public static void mainJob(List<String> strList, int wLen,
			List<Set<Integer>> idSetList, List<Set<Long>> encodeSetList) {

		for (int i = 0; i < strList.size() / 2; i++) {
			Set<Long> encodeSet = stringEncode(strList.get(i * 2), wLen);
			encodeSet.addAll(stringEncode(strList.get(i * 2 + 1), wLen));

			boolean flag = false;

			for (int j = 0; j < encodeSetList.size(); j++) {
				for (Long l : encodeSet) {
					if (encodeSetList.get(j).contains(l)) {
						encodeSetList.get(j).addAll(encodeSet);
						idSetList.get(j).add(i);
						flag = true;
						break;
					}
				}
				if (flag) {
					break;
				}
			}

			if (!flag) {
				Set<Integer> idSet = new TreeSet<Integer>();
				idSet.add(i);
				idSetList.add(idSet);
				encodeSetList.add(encodeSet);
			}
		}

		Collections.sort(idSetList, new Comparator<Set<Integer>>() {
			public int compare(Set<Integer> a1, Set<Integer> a2) {
				return a2.size() - a1.size(); // assumes you want biggest to
												// smallest
			}
		});
	}

	public static List<Set<Integer>> encodeMainJob(List<String> strList, int wLen)
			throws IOException {
		List<Set<Integer>> idSetList = new ArrayList<Set<Integer>>();
		List<Set<Long>> encodeSetList = new ArrayList<Set<Long>>();

		for (int i = 0; i < strList.size() / 2; i++) {
			Set<Long> encodeSet = stringEncode(strList.get(i * 2), wLen);
			encodeSet.addAll(stringEncode(strList.get(i * 2 + 1), wLen));

			boolean flag = false;

			for (int j = 0; j < encodeSetList.size(); j++) {
				for (Long l : encodeSet) {
					if (encodeSetList.get(j).contains(l)) {
						encodeSetList.get(j).addAll(encodeSet);
						idSetList.get(j).add(i);
						flag = true;
						break;
					}
				}
				if (flag) {
					break;
				}
			}

			if (!flag) {
				Set<Integer> idSet = new TreeSet<Integer>();
				idSet.add(i);
				idSetList.add(idSet);
				encodeSetList.add(encodeSet);
			}
		}

		Collections.sort(idSetList, new Comparator<Set<Integer>>() {
			public int compare(Set<Integer> a1, Set<Integer> a2) {
				return a2.size() - a1.size(); // assumes you want biggest to
												// smallest
			}
		});

		return idSetList;
	}
	
	public static void mainJob(List<String> strList, int wLen)
			throws IOException {
		List<Set<Integer>> idSetList = new ArrayList<Set<Integer>>();
		List<Set<Long>> encodeSetList = new ArrayList<Set<Long>>();

		for (int i = 0; i < strList.size() / 2; i++) {
			Set<Long> encodeSet = stringEncode(strList.get(i * 2), wLen);
			encodeSet.addAll(stringEncode(strList.get(i * 2 + 1), wLen));

			boolean flag = false;

			for (int j = 0; j < encodeSetList.size(); j++) {
				for (Long l : encodeSet) {
					if (encodeSetList.get(j).contains(l)) {
						encodeSetList.get(j).addAll(encodeSet);
						idSetList.get(j).add(i);
						flag = true;
						break;
					}
				}
				if (flag) {
					break;
				}
			}

			if (!flag) {
				Set<Integer> idSet = new TreeSet<Integer>();
				idSet.add(i);
				idSetList.add(idSet);
				encodeSetList.add(encodeSet);
			}
		}

		Collections.sort(idSetList, new Comparator<Set<Integer>>() {
			public int compare(Set<Integer> a1, Set<Integer> a2) {
				return a2.size() - a1.size(); // assumes you want biggest to
												// smallest
			}
		});

		for (int i = 0; i < idSetList.size(); i++) {
			for (Integer temp : idSetList.get(i)) {
				System.out.print("\t");
				System.out.print(temp);
			}
			System.out.println();
		}
	}

	public static void decode(Long l) {
		Long mask = 0x3l;
		for (int i = 31; i >= 0; i--) {
			long b = (l >> (i * 2)) & mask;
			char c = 'A';
			if (b == 1) {
				c = 'C';
			} else if (b == 2) {
				c = 'G';
			} else if (b == 3) {
				c = 'T';
			}
			System.out.print(c);
		}
		System.out.println();
	}

	public static void main(String[] args) throws Exception {
		// TODO Auto-generated method stub
		if (args.length < 1) {
			printUsage();
		} else {
			BufferedReader br = new BufferedReader(new FileReader(new File(
					args[0])));
			List<String> strList = new ArrayList<String>();
			String line = null;
			while (null != (line = br.readLine())) {
				strList.add(line);
			}
			br.close();

			mainJob(strList, 32);
		}
	}

}
