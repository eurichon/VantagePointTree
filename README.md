# VantagePointTree
An implementation of the vantage point tree algorithm sequentially and in parallel

Known issues: When getInner(T) or getOuter(T) are called instead of only return the next child object it transform the T to that object 
too. This is because of the array implementation that we chose to model this tree.That problem can be tackled by copying and storing
the T object before any calls of these functions so we can retrive it.
