# yMap

yMap is a python based fast and robust automated method to map large yeast variants to 
            
            proteins post translational modifications 
            
            and other functional regions like:
                            
                        proteins domains, 
                            
                        proteins-DNA binding domains, 
                        
                        proteins structural regions, 
                
                        proteins active and binding sites 
                
                        proteins network visualisation. 


In a user friendly was, it generates "final-report" file to report all the Non-synonymous 
mutations that overlaps or falls inside the above mentioned proteins functional regions. 


    Contents

            scripts: ymap.py

            README: Introduction and details of functions/methods and data with a Manual - with details of how to run yMap 

            Supported data: data/PTMcode + PTMfunc

            Example data: examples mutations files




    Dependencies 
   
    yMap depends on:
            
            python 2.7
            Orange bioinformatics (http://pythonhosted.org/Orange-Bioinformatics/)
    


    Usage
    
            os.chdir(../yMap)
            import ymap
            from ymap import *
            data()
            mutation_types_file()
            ymap()
        
        
    Documentation     
    #For complete instructions see README file

    
#This work in supported by KU Leuven. 
