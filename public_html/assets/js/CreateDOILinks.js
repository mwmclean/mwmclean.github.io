$(document).ready(function(){
    var DOI = document.registerElement('d-o-i', {
      prototype: Object.create(HTMLAnchorElement.prototype),
      extends: 'a'
    });

    $("a").each(function(){
	if ($(this).attr("is") == "d-o-i"){
	  <!-- this.innerHTML = "DOI: " + this.attr("data-doi") + "."; -->
	  this.text = "DOI: " + $(this).attr("data-doi") + ".";
	  $(this).attr("href", "http://dx.doi.org/" + $(this).attr("data-doi"));
	}
    });
});
// function CreateDOILinks() {
//     var DOI = document.registerElement('d-o-i', {
//       prototype: Object.create(HTMLAnchorElement.prototype),
//       extends: 'a'
//     });

//     $("a").each(function(){
// 	if (a.attr("is") == "d-o-i"){
// 	    this.innerHTML = "DOI: " + this.attr("data-doi") + ".";
// 	    this.href = "http://dx.doi.org/" + this.attr("data-doi");
// 	}
//     });
// } 

