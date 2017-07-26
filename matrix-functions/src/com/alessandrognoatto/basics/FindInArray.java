package com.alessandrognoatto.basics;


public class FindInArray
{
	/**
	 * Given an array of double, return the index corresponding to the first element which is less
	 * than value
	 * @param myArray
	 * @param value
	 * @return the first index i such that myArray[i]<value
	 */
	public static int findFirstLessThan(double[] myArray, double value){
		
	    int i = 0;
        boolean found = false;
        
        for ( i = 0; i < myArray.length; i++)
        {
              if (myArray[ i ]  <= value)
              {
                       found = true;     
                       break;
              }
         }
        
        if (found)
        {
        		return i;
        }
        else{
        		System.out.println("No index found");
        		return -1;
        }
        
	}
	
	/**
	 * Given an array of double, return the index corresponding to the first element which is greater
	 * than value
	 * @param myArray
	 * @param value
	 * @return the first index i such that myArray[i]>value
	 */
	public static int findFirstGreaterThan(double[] myArray, double value){
		
	    int i = 0;
        boolean found = false;
        
        for ( i = 0; i < myArray.length; i++)
        {
              if (myArray[ i ]  >= value)
              {
                       found = true;     
                       break;
              }
         }
        
        if (found)
        {
        		return i;
        }
        else{
        		System.out.println("No index found");
        		return -1;
        }
        
	}
	
	 
     public static void main(String[ ] args)
     {
           double[ ] numbers = { 12, 13, 27, 33, 23, 31, 22, 6, 87, 16 };
           double key = 7;
           
           int index = FindInArray.findFirstLessThan(numbers, key);
           int index2 = FindInArray.findFirstGreaterThan(numbers, key);
           
           System.out.println("The first element of the array less than "+key+" is at index "+index);
           System.out.println("The first element of the array greater than "+key+" is at index "+index2);
           
           boolean myBool = false;
           System.out.println(myBool);
           myBool = !myBool;
           System.out.println(myBool);
      }
}