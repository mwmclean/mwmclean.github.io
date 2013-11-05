webpage
=======

My personal webpage

TO DO: 
--------

1) select box for selecting journal (possibly distinguish To Appear from Published)  
2) check document outline and fix if necessary  

Changes to webpage
------------------
### 2013-10-10
1) added html5 semantic elements: changed misc element divs to article, divs with id contentDiv to sections, changed div with id header to header element, added footer element  
2) added collapsible paragraphs to research tab using details tag. h4 tags for project headings changed to summary tags. supported in Opera, Chrome, Safarai. not in Firefox and IE  
3) McLean.css changed accordingly for 1). header id select becomes header element selector, contentDiv selector becomes section selector, summary selector added to h4 rule,   
4) removed \<center\> and \<left\> from research section  
	
### 2013-10-29
1) changed styling of tt and code tags to match GitHub  
2) currently selected navDiv made different colour from other navDivs  
3) In light of 2), removed headings from each div -- this will probably mess up document outline  
3) used padding attribute to contentDiv (section element) to remove reliance on blockquote 
4) changed toggle.js to only take one argument and use jquery 

### 2013-11-05
1) Made navigation divs have rounded border (border-radius attribute)
2) fix issue with underline extending below border of navigation divs (line height and margin attributes)


