$(document).ready(function() {
  //Calls the tocify method on your HTML div.
  $("#toc").tocify();

  //Call the tocify method with options
  $('#toc').tocify({
    showEffect: "fadeIn"
    scrollTo: 50,
    smoothScroll: false
  });
});


$('body').scrollspy({ target: '#toc' })
