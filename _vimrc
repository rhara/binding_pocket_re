set nocompatible

filetype off

set viminfo='50,\"1000,:0,n~/.viminfo

set showmatch      " Show matching brackets.
set ignorecase     " Do case insensitive matching
set smartcase      " Do smart case matching
set incsearch      " Incremental search
set autowrite      " Automatically save before commands like :next and :make
set hidden         " Hide buffers when they are abandoned
set mouse=a        " Enable mouse usage (all modes)

set wildmenu
set wildmode=list:full
set wildignore=*.o,*.obj,*.pyc,*.a,*.so,*.dynlib,*.oeb

set wrapscan

set number
set title
set smartindent
set smarttab
set list
set listchars=tab:»-,trail:-,eol:↲,extends:»,precedes:«,nbsp:%
set nrformats-=octal
set virtualedit=block
set whichwrap=b,s,[,],<,>
set backspace=indent,eol,start
set colorcolumn=100
set tabstop=4
set shiftwidth=4
set expandtab

syntax on

inoremap <silent> jj <ESC>

filetype plugin indent on
