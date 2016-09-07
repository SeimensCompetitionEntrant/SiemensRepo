
/**
 *  Prompt.java
 *  Provides utilities for user input.
 *  "Enhances" the Scanner class, so that
 *  our programs can recover from "bad" input,
 *  and also provide a way to limit numerical
 *  input to a range of values.  Methods for
 *  reading in Strings, ints, and doubles.
 *  @author Your Name Here
 *  @version 1.0
 *  @since 8/20/2015
 */
import java.util.Scanner;
public class Prompt
{
	/**
	 *  Prompts the user and picks up a String.
	 *  @param ask       The String prompt to be displayed to the user.
	 *  @return          The String entered by the user.
	 */
	public static String getString (String ask)
	{
		Scanner keyboard = new Scanner(System.in);
		System.out.print(ask);
		String input = keyboard.nextLine();
		return input;
	}

	/**
	 *  Prompts the user and picks up an int.  Checks for
	 *  "bad" input and reprompts if not an int.
	 *  @param ask       The String prompt to be displayed to the user.
	 *  @return          The int entered by the user.
	 */
	public static int getInt (String ask)
	{
		boolean badinput = false;
		String input = new String("");
		int value = 0;
		do
		{
			badinput = false;
			input = getString(ask);
			try
			{
				value = Integer.parseInt(input);
			}
			catch(NumberFormatException e)
			{
				badinput = true;
			}

		}
		while(badinput);
		return value;
	}

	/**
	 *  Prompts the user and picks up an int.  Checks for
	 *  "bad" input and reprompts if not an int.  Also checks
	 *  for input within a given range, and reprompts if outside
	 *  that range.
	 *  @param ask       The String prompt to be displayed to the user.
	 *  @param min       The minimum integer value to be allowed as input.
	 *  @param max       The maximum integer value to be allowed as input.
	 *  @return          The int entered by the user.
	 */
	public static int getInt (String ask, int min, int max)
	{
		int value = 0;
		do
		{
			value = getInt(ask + " (from " + min + " to " + max + "): ");
		}
		while(value < min || value > max);
		return value;
	}

	public static double getDouble(String ask){
		String input = new String("");
		boolean badinput = false;
		double value = 0.0;

		do{
			badinput = false;
			input = getString(ask);
			try{
				value = Double.parseDouble(input);
			}
			catch(NumberFormatException e){
				badinput = true;
			}
		}
		while(badinput);
		return value;

	}

	public static double getDouble(String ask,double min, double max){
		double value = 0.0;
		do{
			value = getDouble(ask + " from " + min + " to " + max + ".");
		}
		while(value<min || value>max);
		return value;
	}

	public static char getChar(String ask){
		String input = new String("");
		boolean badinput = false;
		char value = ' ';
		do{
			badinput = false;
			input = getString(ask);
			if(input.length()!=1){
				badinput = true;
			}
			else{
				try{
					value = input.charAt(0);
				}
				catch(Exception e){
					System.out.println(e);
				}
			}
		}
		while(badinput);

		return value;

	}

	public static char getChar(String ask,char [] valids){
		char value = ' ';
		boolean badinput = true;
		do{
			value = getChar(ask);
			for(int i = 0; i<valids.length; i++){
				if(value == valids[i]){
					badinput = false;
				}
			}
		}while(badinput);
		return value;
	}
}
