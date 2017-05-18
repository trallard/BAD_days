//layon sidebar
(function(document) {
        var toggle = document.querySelector('.sidebar-toggle');
        var sidebar = document.querySelector('#sidebar');
        var checkbox = document.querySelector('#sidebar-checkbox');
        document.addEventListener('click', function(e) {
          var target = e.target;
          if(!checkbox.checked ||
             sidebar.contains(target) ||
             (target === checkbox || target === toggle)) return;
          checkbox.checked = false;
        }, false);
      })(document);

// table of contents
$('body').scrollspy({ target: '#toc' })


//MD sidebar

// SideNav Options
$('.button-collapse').sideNav({
  edge: 'right', // Choose the horizontal origin
  closeOnClick: true // Closes side-nav on <a> clicks, useful for Angular/Meteor
});



// Show sideNav
$('.button-collapse').sideNav('show');
// Hide sideNav
$('.button-collapse').sideNav('hide');
